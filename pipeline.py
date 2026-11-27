"""QIIME 2 CLI workflow. Commands are argument lists, never shell strings."""
import csv, re, json, shlex, subprocess, os, shutil
from pathlib import Path

RANKS = ['Kingdom','Phylum','Class','Order','Family','Genus','Species']
FASTQ = re.compile(r'\.(?:fastq|fq)(?:\.gz)?$', re.I)
PAIR = re.compile(r'^(.*?)(?:_S\d+)?(?:_L\d{3})?[_\-.]R?([12])(?:_\d{3})?$', re.I)

def manifest(inputs, out, mode, supplied=None):
    files = sorted(p.name for p in inputs.iterdir() if p.is_file() and FASTQ.search(p.name))
    paired = mode == 'paired'
    columns = ['sample-id','forward-absolute-filepath','reverse-absolute-filepath'] if paired else ['sample-id','absolute-filepath']
    rows = {}
    if supplied:
        with (inputs / supplied).open(encoding='utf-8-sig', newline='') as f:
            reader=csv.DictReader(f, delimiter='\t'); required=['sample-id','forward','reverse'] if paired else ['sample-id','file']
            if not set(required)<=set(reader.fieldnames or []): raise ValueError('Manifest TSV needs '+', '.join(required))
            for r in reader:
                s=r['sample-id'].strip(); names=[r[k].strip() for k in required[1:]]
                if s in rows: raise ValueError('Duplicate sample ID: '+s)
                if any(n not in files for n in names): raise ValueError('Manifest references missing FASTQ: '+str(names))
                rows[s]=names
    else:
        for name in files:
            stem=FASTQ.sub('',name)
            if paired:
                m=PAIR.fullmatch(stem)
                if not m: raise ValueError('Cannot identify R1/R2 for '+name+'; provide manifest TSV.')
                s, read=m.groups(); row=rows.setdefault(s,[None,None]); idx=int(read)-1
                if row[idx]: raise ValueError('Multiple lanes/reads for '+s+'; combine lanes per sample before upload.')
                row[idx]=name
            else:
                if stem in rows: raise ValueError('Duplicate sample '+stem)
                rows[stem]=[name]
    if not rows: raise ValueError('No FASTQ samples found.')
    seen=set()
    with out.open('w',newline='') as f:
        writer=csv.writer(f,delimiter='\t'); writer.writerow(columns)
        for s,names in sorted(rows.items()):
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',s) or any(n is None for n in names): raise ValueError('Invalid ID or missing mate: '+s)
            if len(set(names))!=len(names) or seen.intersection(names): raise ValueError('FASTQ reused in manifest: '+s)
            seen.update(names); writer.writerow([s]+[str((inputs/n).resolve()) for n in names])
    return sorted(rows)

def metadata_file(inputs, name, out, ids=None):
    if not name:
        if ids is None: raise ValueError('Metadata required for artifact workflow.')
        with out.open('w') as f:
            f.write('sample-id\tSample\n'); f.writelines(f'{s}\t{s}\n' for s in ids)
        return
    import pandas as pd
    p=inputs/name; df=pd.read_csv(p,sep=',' if p.suffix.lower()=='.csv' else '\t',dtype=str,keep_default_na=False)
    if df.shape[1]<2: raise ValueError('Metadata needs sample ID and at least one variable.')
    first=df.columns[0]; actual=df[~df[first].str.startswith('#q2:')][first]
    if actual.duplicated().any() or (actual.str.strip()=='').any(): raise ValueError('Metadata contains duplicate or empty IDs.')
    if ids and not set(ids)<=set(actual): raise ValueError('Metadata missing sample IDs: '+', '.join(sorted(set(ids)-set(actual))[:15]))
    df=df.rename(columns={first:'sample-id'}); df.to_csv(out,sep='\t',index=False)

def execute(project, run, cfg, report):
    if cfg['mode']=='ont':
        from ont import execute_ont
        return execute_ont(project,run,cfg,report)
    inp=project/'inputs'; log=run/'commands.log'; n=min(int(cfg.get('threads',4)),int(os.environ.get('MAX_THREADS','16')))
    def command(label,*args):
        cmd=['qiime',*map(str,args)]; report(label)
        with log.open('a') as f:
            f.write('\n$ '+shlex.join(cmd)+'\n'); f.flush()
            subprocess.run(cmd, cwd=run, stdout=f, stderr=subprocess.STDOUT, check=True)
    def q(name):
        if not name or Path(name).name!=name or not (inp/name).is_file(): raise ValueError('Missing input: '+str(name))
        return inp/name
    mode=cfg['mode']; stage=cfg.get('stage','full'); md=run/'metadata.tsv'; table=run/'table.qza'; seqs=run/'rep-seqs.qza'; tax=run/'taxonomy.qza'; tree=run/'rooted-tree.qza'
    if mode in ('single','paired','ccs'):
        ids=manifest(inp,run/'manifest.tsv',mode,cfg.get('manifest'))
        demux=run/'demux.qza'; paired=mode=='paired'
        command('Import FASTQ','tools','import','--type','SampleData[PairedEndSequencesWithQuality]' if paired else 'SampleData[SequencesWithQuality]','--input-format','PairedEndFastqManifestPhred33V2' if paired else 'SingleEndFastqManifestPhred33V2','--input-path',run/'manifest.tsv','--output-path',demux)
        command('Read quality report','demux','summarize','--i-data',demux,'--o-visualization',run/'quality.qzv')
        if stage=='quality': return
        f=cfg.get('primer_f',''); r=cfg.get('primer_r','')
        if mode!='ccs' and (f or r):
            if not f or (paired and not r) or (not paired and r): raise ValueError('Supply forward primer and reverse primer only for paired reads.')
            trimmed=run/'trimmed.qza'; opts=['--p-front-f',f,'--p-front-r',r] if paired else ['--p-front',f]
            command('Remove primers','cutadapt','trim-paired' if paired else 'trim-single','--i-demultiplexed-sequences',demux,*opts,'--p-cores',n,'--o-trimmed-sequences',trimmed); demux=trimmed
        if mode=='paired':
            method='denoise-paired'; opts=['--p-trunc-len-f',cfg['trunc_f'],'--p-trunc-len-r',cfg['trunc_r'],'--p-trim-left-f',cfg['trim_f'],'--p-trim-left-r',cfg['trim_r']]
        elif mode=='single':
            method='denoise-single'; opts=['--p-trunc-len',cfg['trunc_f'],'--p-trim-left',cfg['trim_f']]
        else:
            if not f or not r: raise ValueError('CCS requires both primer sequences.')
            method='denoise-ccs'; opts=['--p-front',f,'--p-adapter',r,'--p-min-len',cfg['min_len'],'--p-max-len',cfg['max_len']]
        command('DADA2 denoising and chimera removal','dada2',method,'--i-demultiplexed-seqs',demux,*opts,'--p-n-threads',n,'--o-table',table,'--o-representative-sequences',seqs,'--o-denoising-stats',run/'denoising-stats.qza','--o-base-transition-stats',run/'base-transition-stats.qza')
        command('DADA2 error model plot','dada2','plot-base-transitions','--i-base-transition-stats',run/'base-transition-stats.qza','--o-visualization',run/'base-transition-stats.qzv')
        command('Denoising statistics','metadata','tabulate','--m-input-file',run/'denoising-stats.qza','--o-visualization',run/'denoising-stats.qzv')
        metadata_file(inp,cfg.get('metadata'),md,ids)
    else:
        tables=[q(x) for x in cfg.get('tables',[])]; imported=[]
        for i,t in enumerate(tables):
            if t.suffix.lower()=='.biom':
                target=run/f'imported-table-{i+1}.qza'; command('Import BIOM frequency table','tools','import','--input-path',t,'--type','FeatureTable[Frequency]','--output-path',target); imported.append(target)
            else: imported.append(t)
        tables=imported; sequences=[q(x) for x in cfg.get('sequences',[])]; taxonomies=[q(x) for x in cfg.get('taxonomies',[])]
        if not tables: raise ValueError('Choose at least one feature table.')
        if len(tables)>1: command('Merge feature tables','feature-table','merge','--i-tables',*tables,'--p-overlap-method','sum' if cfg.get('sum_overlap') else 'error_on_overlapping_sample','--o-merged-table',table)
        else: shutil.copyfile(tables[0],table)
        if sequences:
            if len(sequences)>1: command('Merge representative sequences','feature-table','merge-seqs','--i-data',*sequences,'--o-merged-data',seqs)
            else: shutil.copyfile(sequences[0],seqs)
        if taxonomies:
            if len(taxonomies)>1: command('Merge taxonomy','feature-table','merge-taxa','--i-data',*taxonomies,'--o-merged-data',tax)
            else: shutil.copyfile(taxonomies[0],tax)
        if cfg.get('tree'): shutil.copyfile(q(cfg['tree']),tree)
        metadata_file(inp,cfg.get('metadata'),md)
    command('Feature table summary','feature-table','summarize','--i-table',table,'--m-sample-metadata-file',md,'--o-visualization',run/'table.qzv')
    if cfg.get('classifier'):
        if not seqs.exists(): raise ValueError('Classifier requires representative sequences.')
        if tax.exists(): tax.unlink()
        command('Taxonomy assignment','feature-classifier','classify-sklearn','--i-classifier',q(cfg['classifier']),'--i-reads',seqs,'--p-n-jobs',n,'--o-classification',tax)
    if tax.exists():
        if cfg.get('exclude_organelle'):
            filtered=run/'filtered-table.qza'; command('Filter mitochondria/chloroplast','taxa','filter-table','--i-table',table,'--i-taxonomy',tax,'--p-exclude','mitochondria,chloroplast','--o-filtered-table',filtered); table=filtered
        command('Taxonomic composition','taxa','barplot','--i-table',table,'--i-taxonomy',tax,'--m-metadata-file',md,'--o-visualization',run/'taxa-barplot.qzv')
        for level in cfg.get('levels',[]):
            if level not in range(1,8): raise ValueError('Invalid taxonomic rank')
            target=run/f'level-{level}.qza'; command('Collapse '+RANKS[level-1],'taxa','collapse','--i-table',table,'--i-taxonomy',tax,'--p-level',level,'--o-collapsed-table',target)
            export_table(command,target,run/f'level-{level}')
    export_table(command,table,run/'feature-table')
    if seqs.exists(): command('Export representative sequences','tools','export','--input-path',seqs,'--output-path',run/'sequences')
    if tax.exists(): command('Export taxonomy','tools','export','--input-path',tax,'--output-path',run/'taxonomy')
    if cfg.get('phylogeny') and not tree.exists():
        if not seqs.exists(): raise ValueError('Phylogeny requires sequences or uploaded tree.')
        command('MAFFT alignment and FastTree','phylogeny','align-to-tree-mafft-fasttree','--i-sequences',seqs,'--p-n-threads',n,'--o-alignment',run/'aligned.qza','--o-masked-alignment',run/'masked.qza','--o-tree',run/'unrooted-tree.qza','--o-rooted-tree',tree)
    if cfg.get('diversity'):
        depth=int(cfg.get('depth',0))
        if depth<=0: raise ValueError('Choose a positive rarefaction depth from feature table summary.')
        opts=['--i-phylogeny',tree] if tree.exists() else []
        command('Alpha rarefaction','diversity','alpha-rarefaction','--i-table',table,*opts,'--p-max-depth',depth,'--m-metadata-file',md,'--o-visualization',run/'alpha-rarefaction.qzv')
        command('Core alpha/beta diversity and PCoA','diversity','core-metrics-phylogenetic' if tree.exists() else 'core-metrics','--i-table',table,*opts,'--p-sampling-depth',depth,'--m-metadata-file',md,'--output-dir',run/'diversity')
        group=cfg.get('group','').strip()
        if group:
            for v in sorted((run/'diversity').glob('*_vector.qza')):
                command('Alpha group significance','diversity','alpha-group-significance','--i-alpha-diversity',v,'--m-metadata-file',md,'--o-visualization',v.with_suffix('.qzv'))
            for v in sorted((run/'diversity').glob('*_distance_matrix.qza')):
                for method in ('permanova','permdisp'):
                    command(method.upper(),'diversity','beta-group-significance','--i-distance-matrix',v,'--m-metadata-file',md,'--m-metadata-column',group,'--p-method',method,'--o-visualization',v.with_name(v.stem+'-'+method+'.qzv'))

def export_table(command, artifact, folder):
    command('Export '+artifact.stem,'tools','export','--input-path',artifact,'--output-path',folder)
    with (folder.parent/'commands.log').open('a') as f:
        subprocess.run(['biom','convert','-i',str(folder/'feature-table.biom'),'-o',str(folder/'counts.tsv'),'--to-tsv'],stdout=f,stderr=subprocess.STDOUT,check=True)
