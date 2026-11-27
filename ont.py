"""ONT full-length 16S workflow, with streaming QC and Emu taxon profiling."""
import csv, gzip, math, os, re, shlex, subprocess, shutil
from pathlib import Path
from collections import defaultdict
import pandas as pd
import numpy as np

FASTQ=re.compile(r'\.(?:fastq|fq)(?:\.gz)?$', re.I)
RANKS=['Kingdom','Phylum','Class','Order','Family','Genus','Species']
PREFIX=['k','p','c','o','f','g','s']
ERROR=[10**(-i/10) for i in range(94)]

def sample_map(inp, manifest=None):
    files={p.name:p for p in inp.iterdir() if p.is_file() and FASTQ.search(p.name)}
    rows=defaultdict(list); used=set()
    if manifest:
        with (inp/manifest).open(encoding='utf-8-sig',newline='') as f:
            reader=csv.DictReader(f,delimiter='\t')
            if not {'sample-id','file'}<=set(reader.fieldnames or []): raise ValueError('ONT manifest needs sample-id and file columns')
            for row in reader:
                s=row['sample-id'].strip(); name=row['file'].strip()
                if name not in files or name in used: raise ValueError('Missing or repeated FASTQ: '+name)
                rows[s].append(files[name]); used.add(name)
    else:
        for name in sorted(files):
            s=FASTQ.sub('',name)
            if s in rows: raise ValueError('Duplicate sample filename stem; use a manifest.')
            rows[s].append(files[name])
    if not rows: raise ValueError('No demultiplexed FASTQ files selected.')
    for s in rows:
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',s): raise ValueError('Invalid sample ID: '+s)
    return dict(sorted(rows.items()))

def records(paths):
    for p in paths:
        opener=gzip.open if p.name.lower().endswith('.gz') else open
        with opener(p,'rb') as f:
            while True:
                h=f.readline()
                if not h: break
                s=f.readline().rstrip(b'\r\n'); plus=f.readline(); q=f.readline().rstrip(b'\r\n')
                if not h.startswith(b'@') or not plus.startswith(b'+') or not s or len(s)!=len(q) or min(q)<33 or max(q)>126:
                    raise ValueError('Malformed four-line Phred+33 FASTQ: '+p.name)
                yield h.rstrip(b'\r\n'),s,plus.rstrip(b'\r\n'),q

def scan(paths,out=None,min_len=0,max_len=10**9,min_q=0):
    n=bases=kept=kept_bases=length_fail=quality_fail=0; qtotal=0.0
    hist=defaultdict(int)
    handle=gzip.open(out,'wb',compresslevel=3) if out else None
    try:
        for h,s,p,q in records(paths):
            l=len(s); score=-10*math.log10(sum(ERROR[x-33] for x in q)/l)
            n+=1; bases+=l; qtotal+=score; hist[(l//100)*100]+=1
            if not min_len<=l<=max_len: length_fail+=1; continue
            if score<min_q: quality_fail+=1; continue
            kept+=1; kept_bases+=l
            if handle: handle.write(b'\n'.join((h,s,p,q))+b'\n')
    finally:
        if handle: handle.close()
    return dict(reads=n,bases=bases,mean_length=bases/n if n else 0,mean_read_q=qtotal/n if n else 0,retained_reads=kept,retained_bases=kept_bases,length_fail=length_fail,quality_fail=quality_fail),hist

def profile_tables(profiles):
    counts={}; rel={}; lineages={}; excluded=[]
    for sample,path in profiles.items():
        d=pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False)
        d.columns=[c.strip().lower() for c in d.columns]
        required={'tax_id','abundance','estimated counts'}
        if not required<=set(d.columns): raise ValueError('Unexpected Emu output columns: '+str(d.columns.tolist()))
        c={}; a={}
        for _,r in d.iterrows():
            tid=r['tax_id'].strip()
            if not tid.isdigit() or int(tid)==0:
                excluded.append(dict(sample=sample,tax_id=tid,abundance=r['abundance'],estimated_counts=r['estimated counts'])); continue
            if tid in c: raise ValueError('Duplicate tax ID: '+tid)
            count=float(r['estimated counts']); abundance=float(r['abundance'])
            if not math.isfinite(count) or not math.isfinite(abundance) or count<0 or abundance<0: raise ValueError('Invalid Emu abundance')
            if count==0 and abundance==0: continue
            vals=[]
            for rank,prefix in zip(RANKS,PREFIX):
                val=r.get(rank.lower(),'')
                if rank=='Kingdom' and val in ('','-','NA','nan'): val=r.get('superkingdom',r.get('domain',''))
                val='' if val in ('','-','NA','nan','None') else val.replace(';',',').strip()
                vals.append(prefix+'__'+val)
            lineage='; '.join(vals)
            if tid in lineages and lineages[tid]!=lineage: raise ValueError('Inconsistent taxonomy for '+tid)
            lineages[tid]=lineage;c[tid]=count;a[tid]=abundance
        counts[sample]=pd.Series(c,dtype=float);rel[sample]=pd.Series(a,dtype=float)
    c=pd.DataFrame(counts).fillna(0); a=pd.DataFrame(rel).fillna(0)
    c.index.name=a.index.name='Feature ID'; a=a.div(a.sum(axis=0).replace(0,np.nan),axis=1).fillna(0)
    taxonomy=pd.Series(lineages,name='Taxon').rename_axis('Feature ID')
    return c,a,taxonomy,pd.DataFrame(excluded,columns=['sample','tax_id','abundance','estimated_counts'])

def execute_ont(project,run,cfg,report):
    from pipeline import metadata_file
    from tables import collapse
    inp=project/'inputs'; samples=sample_map(inp,cfg.get('manifest')); log=run/'commands.log'
    n=min(int(cfg.get('threads',4)),int(os.environ.get('MAX_THREADS','16')))
    def cmd(label,args):
        report(label); args=list(map(str,args))
        with log.open('a') as f:
            f.write('\n$ '+shlex.join(args)+'\n'); f.flush(); subprocess.run(args,cwd=run,stdout=f,stderr=subprocess.STDOUT,check=True)
    db=Path(os.environ.get('EMU_DATABASE_DIR','/data/emu_db'))
    full=cfg['stage']=='full'; f=cfg.get('primer_f',''); r=cfg.get('primer_r','')
    if bool(f)!=bool(r): raise ValueError('Provide both primers or neither.')
    if full and not all((db/x).is_file() for x in ('species_taxid.fasta','taxonomy.tsv')): raise ValueError('Install the Emu database on the worker: '+str(db))
    minlen=int(cfg.get('min_len',1200)); maxlen=int(cfg.get('max_len',1800)); minq=float(cfg.get('min_q',10))
    if not 1<=minlen<=maxlen or not 0<=minq<=60: raise ValueError('Invalid ONT filtering parameters')
    scratch=run/'scratch';scratch.mkdir(); result=run/'emu';result.mkdir(); rows=[]; histrows=[]; profiles={}
    pd.DataFrame([(s,p.name) for s,paths in samples.items() for p in paths],columns=['sample-id','file']).to_csv(run/'ont-manifest.tsv',sep='\t',index=False)
    for s,paths in samples.items():
        report('ONT raw QC: '+s); raw,hist=scan(paths)
        histrows.extend({'sample-id':s,'length_bin_start':k,'reads':v} for k,v in sorted(hist.items()))
        row={'sample-id':s,**{'raw_'+k:v for k,v in raw.items() if k in ('reads','bases','mean_length','mean_read_q')}}
        source=paths; merged=scratch/'merged.fastq.gz'; trimmed=scratch/'trimmed.fastq.gz'; filtered=scratch/'filtered.fastq.gz'
        try:
            if f:
                report('Combine chunks: '+s); scan(paths,merged)
                rc=r.upper().translate(str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN'))[::-1]
                cmd('Trim both terminal primers: '+s,[os.environ.get('CUTADAPT_BIN','/opt/emu/bin/cutadapt'),'--revcomp','--discard-untrimmed','-j',n,'-g','^'+f.upper()+'...'+rc+'$','-o',trimmed,merged]); source=[trimmed]
            report('Length and quality filtering: '+s); stats,_=scan(source,filtered if full else None,minlen,maxlen,minq)
            row.update(stats);row['primer_retained_reads']=stats['reads']; rows.append(row)
            pd.DataFrame(rows).to_csv(run/'ont-qc.tsv',sep='\t',index=False)
            pd.DataFrame(histrows).to_csv(run/'read-length-histogram.tsv',sep='\t',index=False)
            if full and stats['retained_reads']:
                cmd('Emu taxonomic profiling: '+s,[os.environ.get('EMU_BIN','/opt/emu/bin/emu'),'abundance',filtered,'--type','map-ont','--db',db,'--threads',n,'--keep-counts','--output-dir',result,'--output-basename',s]); path=result/(s+'_rel-abundance.tsv')
                if not path.is_file(): raise ValueError('Emu did not create expected profile: '+s)
                profiles[s]=path
        finally:
            for p in (merged,trimmed,filtered): p.unlink(missing_ok=True)
    shutil.rmtree(scratch)
    if not full:return
    if not profiles:raise ValueError('No sample retained reads after filtering; see QC.')
    metadata_file(inp,cfg.get('metadata'),run/'metadata.tsv',list(samples))
    c,a,tax,excluded=profile_tables(profiles)
    excluded.to_csv(run/'unclassified-profile-rows.tsv',sep='\t',index=False)
    c=c.reindex(columns=list(samples),fill_value=0);a=a.reindex(columns=list(samples),fill_value=0)
    c.to_csv(run/'estimated-counts-all-samples.tsv',sep='\t'); a.to_csv(run/'relative-abundance-all-samples.tsv',sep='\t');tax.to_frame().to_csv(run/'taxonomy.tsv',sep='\t')
    keep=a.sum(axis=0)>0
    pd.DataFrame({'sample-id':a.columns,'included_in_composition':keep.values}).to_csv(run/'sample-inclusion.tsv',sep='\t',index=False)
    if not keep.any():raise ValueError('No classified taxa; see Emu outputs and QC.')
    a=a.loc[:,keep];c=c.loc[:,keep]
    for level in cfg.get('levels',range(1,8)):
        for name,table in [('estimated-counts',c),('relative-abundance',a)]: collapse(table,tax,level).to_csv(run/f'{RANKS[level-1].lower()}-{name}.tsv',sep='\t')
    report('Generate descriptive diversity and QIIME 2 relative-frequency artifact')
    p=a.where(a>0);pd.DataFrame({'observed_taxa':(a>0).sum(),'shannon_taxon':-(p*np.log(p)).sum()}).rename_axis('sample-id').to_csv(run/'taxon-alpha-diversity.tsv',sep='\t')
    if a.shape[1]>=2:
        from scipy.spatial.distance import pdist,squareform
        d=squareform(pdist(a.T,metric='braycurtis'));pd.DataFrame(d,index=a.columns,columns=a.columns).to_csv(run/'bray-curtis.tsv',sep='\t')
        J=np.eye(len(d))-np.ones_like(d)/len(d); vals,vecs=np.linalg.eigh(-0.5*J@(d*d)@J);order=np.argsort(vals)[::-1];vals=vals[order];vecs=vecs[:,order];coords=vecs[:,:2]*np.sqrt(np.maximum(vals[:2],0));pd.DataFrame(coords,index=a.columns,columns=['PC1','PC2']).rename_axis('sample-id').to_csv(run/'bray-curtis-pcoa.tsv',sep='\t')
    from biom import Table
    from qiime2 import Artifact
    Artifact.import_data('FeatureTable[RelativeFrequency]',Table(a.to_numpy(),a.index.tolist(),a.columns.tolist())).save(str(run/'ont-relative-frequency.qza'))
    cmd('QIIME 2 taxonomy import',['qiime','tools','import','--type','FeatureData[Taxonomy]','--input-path',run/'taxonomy.tsv','--output-path',run/'ontology-taxonomy.qza'])
    cmd('QIIME 2 taxonomic barplot',['qiime','taxa','barplot2','--i-table',run/'ont-relative-frequency.qza','--i-taxonomy',run/'ontology-taxonomy.qza','--m-metadata-file',run/'metadata.tsv','--o-visualization',run/'ont-taxa-barplot.qzv'])
