# Microbiome Workbench | QIIME 2 + Streamlit

An editable private research web portal, with a Streamlit interface and a separate persistent analysis worker. It is not yet hosted or production-validated. The bundled CSV examples are synthetic, for testing only.

## Architecture and file limits

Browser → Streamlit (settings, plots, status) → FastAPI worker (job queue) → QIIME 2 2026.10 workshop container. The browser uploads FASTQ/QZA directly to the worker in 16 MiB chunks; interrupted uploads can resume by selecting the exact same file. Streamlit's own file uploader is configured to 10,240 MB for table tools, but is memory-bound and should only be used for small tables. The direct uploader defaults to 100 GiB per file and 500 GiB completed/partial inputs per project; change `MAX_FILE_GIB` and `PROJECT_GIB` on the worker. Output files and scratch are additional; there is no automatic output quota or expiry. Monitor disk and clean up projects. Uploading 100 GiB has not been load-tested. File count and processing capacity depend on server resources; there is no claim of unlimited uploads.

## 1. Run both services locally or on an institutional VM

Requirements: Linux x86-64, Docker Engine + Compose, sufficient RAM, storage, internet to pull container and dependencies. Starting hardware suggestion for moderate batches: 16 CPU threads, 64 GB RAM, 1–2 TB SSD; large or many full-length samples may need more. Estimate scratch/output space before each run.

```bash
cp .env.example .env
mkdir -p data/tmp data/projects
# Edit .env: replace both secrets with independent random values (openssl rand -hex 32).
# Set WORKER_PUBLIC_URL to your actual HTTPS worker endpoint for remote browsers.
docker compose build
docker compose up -d
docker compose exec worker qiime info
```

Open http://localhost:8501 in the browser on the machine running Docker (or use SSH port forwarding for both ports). Worker uploader: localhost:8000. Both ports are loopback-bound by default. Set writable permissions on the mounted data directory for the worker container's user if necessary. The QIIME image and runtime are not available in the generation environment: validate `qiime info`, DADA2, classifier and a small end-to-end dataset before using study data. Run `docker compose logs worker` for infrastructure troubleshooting.

## 2. Deploy the Streamlit frontend on streamlit.app

1. Create a private GitHub repository and upload this directory's **contents**, including `.streamlit/config.toml`. Do not commit `.env`, secrets, raw data, or classifiers. GitHub is a code host, not a sequencing data host.
2. Start the worker on a Linux VM with `docker compose up -d --build worker`. Make the worker reachable at an HTTPS domain via a reverse proxy; see example below. Ensure the VM has persistent storage and the worker port is not exposed directly.
3. In [Streamlit Community Cloud](https://share.streamlit.io), create an app from the GitHub repository, select `app.py` as the main file, and add the contents of `.streamlit/secrets.toml.example` to Advanced settings → Secrets, replacing URLs and both secrets with your real values. The `WORKER_TOKEN` must match the worker environment. Restrict app access to authorised users as well as the built-in shared password.
4. Deploy. Sign in, create a project, then click **Open large-file uploader**. FASTQ files will bypass Streamlit entirely. When finished, return to the portal and refresh file list. Choose quality check or automatic analysis.

For example, a Caddy reverse proxy on the VM (DNS and TLS must be configured):

```
qiime-worker.your-domain.example {
    reverse_proxy 127.0.0.1:8000
}
```

Do not enable logging of full query strings: project access keys are capability URLs and grant access to data. HTTPS, private app access, approved institutional hosting, network controls, disk encryption, backups and a data-retention policy are essential before research use. The starter uses a shared password and random project access keys, not institutional SSO, per-user roles, audit trails or multi-tenant isolation. Do not allow untrusted people to upload classifier artifacts; serialized classifiers can execute code when loaded. For strong isolation, use separate trusted-user workers with least privilege. No third-party upload service is required. Do not send identifiable human data to unauthorised hosting services.

Community Cloud is appropriate for the frontend and small table utilities, not for running raw sequencing pipelines within its resource limits. If your institution requires data to remain inside its network, host both services privately instead and use VPN/SSO gateway. A public cloud frontend can see selected inputs/results and metadata filenames; assess governance before use.

## 3. Workflow

- Upload demultiplexed Phred+33 ONT full-length 16S, Illumina single/paired reads, PacBio CCS/HiFi reads, or existing QZA and BIOM v2.1 count tables. Choose the correct mode. For non-ONT modes, multiple sequencing lanes for one sample must be combined upstream; different runs should be denoised separately before merging. Provide a filename manifest for unusual names.
- Run read-quality report first; download QZV and view in QIIME 2 View (use approved local viewing for sensitive studies). Confirm primers and truncation from QC; for paired reads preserve sufficient overlap. For CCS confirm both adapter orientations and length range. For Oxford Nanopore full-length 16S, use the dedicated ONT mode described below; shotgun reads are not supported. Preprocessed count tables can be analysed independently of technology.
- Select full pipeline, optional cutadapt, DADA2 denoising, taxonomy classifier, optional organelle exclusion, taxonomic levels, phylogenetic tree and diversity. No reference database/classifier is bundled. Use a trusted, compatible classifier for your sequencing region and QIIME release. Taxonomy outputs are omitted if no classification or taxonomy input is supplied.
- For a scientifically safer two-stage run, select denoising/taxonomy without diversity first. Review feature table summary, select a rarefaction depth, and reuse `table.qza` (or `filtered-table.qza` if organelles excluded), `rep-seqs.qza`, `taxonomy.qza`, `metadata.tsv`, and rooted tree if available via Results → Reuse. Select these in existing-artifact mode and run diversity. This avoids re-upload and repeats neither read import nor denoising.
- Results: QZA/QZV, BIOM, count TSVs, sequence FASTA, taxonomy TSV, command log and exact run parameters. QIIME core diversity calculates Shannon, observed features, Pielou evenness, Bray–Curtis/Jaccard, and, with a valid rooted tree, Faith PD and weighted/unweighted UniFrac; PCoA/Emperor plots. Group significance is exploratory; repeated measurements, batch effects, covariates and multiple comparisons require study-specific models. Alpha significance explores categorical metadata; beta tests use the specified column. Species assignment is not necessarily reliable even with full-length 16S.
- Independent table utilities: merge by ID; taxonomic collapse from Kingdom to Species; collapse each then merge; sum/mean/median by metadata groups; plot count and relative abundance tables. Keep sample IDs unique; only sum technical replicates. Table operations require data small enough for frontend RAM. Use worker QZA/BIOM workflow for large artifacts. Table utilities accept a taxonomy file with `Feature ID` first column and `Taxon` semicolon-delimited ordered lineage column.
- Use synthetic files in `examples/` for a first frontend test. Merge the two count files; use taxonomy.tsv and Genus rank; optionally aggregate by Group from metadata.tsv. Visualisation accepts the result or a CSV/TSV; download interactive HTML. QIIME native plots download as QZV.

## 4. Validation status and limitations

Python modules compiled and three automated tests passed: filename pairing and missing mates; table merge/overlap, collapse, metadata aggregation and relative abundance; worker authentication, chunk upload/resume, duplicate offset rejection, file listing and deletion. Actual Docker/QIIME execution, Streamlit UI browser rendering, public deployment and 100 GiB stress testing remain unverified. This is a working-source prototype, not a certified production application. The single worker serializes runs; jobs survive browser disconnection but not a worker restart. A failed job remains available for inspection and can be rerun. No differential-abundance, contamination inference, functional prediction, longitudinal modelling or automated primer/truncation optimisation is bundled. QIIME action availability and command flags are based on 2026.10 documentation.

Reference documentation: https://docs.qiime2.org/en/2026.10/ ; https://library.qiime2.org/quickstart/qiime2 ; https://docs.streamlit.io/knowledge-base/deploy/increase-file-uploader-limit-streamlit-cloud ; https://docs.streamlit.io/deploy/streamlit-community-cloud/manage-your-app

## Oxford Nanopore full-length 16S (default mode)

This edition adds a dedicated ONT workflow: basecalled, already demultiplexed Phred+33 FASTQ/FASTQ.GZ → per-sample read QC → optional terminal primer removal in both orientations (Cutadapt) → length and mean-read Phred filtering → Emu 3.6.2 species-level profiling → estimated-count and relative-abundance tables → collapse at all selected ranks → descriptive taxon richness, Shannon and Bray–Curtis PCoA → QIIME 2 relative-frequency artifact and taxonomic plot. No standard Illumina DADA2 denoising is applied to ONT reads. No true ASV inference or standalone chimera detection is implemented. Species assignments are database-dependent; check positive/negative controls and the protocol. This is not a shotgun metagenomics pipeline. Basecalling and barcode demultiplexing must be completed upstream.

### One-time database setup

Download a compatible Emu database using the official instructions at https://github.com/treangenlab/emu#installation (default database archive on https://osf.io/56uf7/). For the official default database, install osfclient locally, then run:

```bash
mkdir -p data/emu_db
cd data/emu_db
osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar
tar -xf emu.tar
# If the archive creates a nested directory, move the two required files to this directory.
ls species_taxid.fasta taxonomy.tsv
cd ../..
docker compose build worker
docker compose up -d
docker compose exec worker /opt/emu/bin/emu --help
docker compose exec worker /opt/emu/bin/cutadapt --version
```

A custom compatible database can be placed in the same directory. Record database version, provenance and file checksums in your study records. The default Emu database is based on a historical NCBI/RRNDB release; alternative SILVA databases have different coverage and taxonomy. Reference data are not bundled. Docker build and actual QIIME/Emu execution require validation on the deployment server.

### Upload and run

1. Upload ONT FASTQ files, with unique filenames. For one FASTQ per sample, sample ID is the filename stem. For multiple chunks per barcode/sample, upload an explicit TSV manifest containing `sample-id` and `file`; sample IDs may repeat, filenames may not. Rename duplicate filenames before upload. Never combine barcodes based on guessing. The manifest may select a subset of project FASTQs. Use `examples/ont-manifest.tsv` as a template.
2. Select Oxford Nanopore full-length 16S in Run analysis. Select the manifest if needed and optional sample metadata. The defaults of 1,200–1,800 bp and mean Q10 are adjustable starting values; check the read distribution, library kit, basecaller and controls. Use quality-only mode first if uncertain. Leave primer removal off if the vendor already trimmed reads; otherwise provide both actual primers in 5′→3′ orientation. The linked primer filter requires both ends and can remove partial reads. The length thresholds apply after trimming.
3. Start analysis. Check Results & logs. QC includes raw and retained read totals, mean lengths/Q and a read-length histogram TSV. Empty or unclassified samples are documented in sample-inclusion.tsv and excluded from compositional diversity. Per-sample Emu output files are retained; unclassified summary rows are exported separately. The relative table is normalized among classified taxa only. Estimated counts are fractional likelihood-based values and are preserved, not rounded or converted to QIIME frequency artifacts. QIIME artifact is `FeatureTable[RelativeFrequency]`, not compatible with count-only artifact mode, count rarefaction or UniFrac. Taxon profiles do not have ASV representative sequences or a phylogenetic tree. Diversity is descriptive, not ASV-level or inferential. Do not sum relative abundances across technical replicates; use compatible estimated-count tables or combine raw chunks before profiling.
4. Download the taxon abundance tables (Kingdom through Species), estimated-counts-all-samples.tsv, relative-abundance-all-samples.tsv, taxonomy.tsv, and QZA/QZV. Use the independent table tools and visualisation tab with these TSVs. Columns for samples without classified taxa are zero in the all-samples exports; treat them as failed/excluded, not biological zeros. Use sample-inclusion.tsv. Raw observed read totals are in ont-qc.tsv, not Emu's estimated-count table.

Emu reference: Curry et al. Nature Methods (2022), https://doi.org/10.1038/s41592-022-01520-4. ONT protocol: https://doi.org/10.1002/cpz1.978. QIIME relative-frequency import: https://docs.qiime2.org/en/stable/how-to-guides/how-to-import/.
