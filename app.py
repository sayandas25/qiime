import os, io, hmac, json
from urllib.parse import quote, urlencode
import requests
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from tables import load_counts,merge_counts,collapse,aggregate,relative

st.set_page_config(page_title='Microbiome Workbench',page_icon='🧬',layout='wide')
RANKS=['Kingdom','Phylum','Class','Order','Family','Genus','Species']
def setting(name,default=''):
    if name in os.environ:return os.environ[name]
    try:return st.secrets.get(name,default)
    except (FileNotFoundError,KeyError):return default
BASE=setting('WORKER_URL','').rstrip('/'); PUBLIC=setting('WORKER_PUBLIC_URL',BASE).rstrip('/'); TOKEN=setting('WORKER_TOKEN'); PASSWORD=setting('PORTAL_PASSWORD')
if not PASSWORD:
    st.error('Set PORTAL_PASSWORD in Streamlit secrets or environment before using this private portal.'); st.stop()
if not st.session_state.get('logged_in'):
    st.title('Microbiome Workbench'); st.caption('Private research workspace · QIIME 2')
    with st.form('login'):
        pwd=st.text_input('Workspace password',type='password'); submit=st.form_submit_button('Sign in')
    if submit:
        if hmac.compare_digest(pwd,PASSWORD): st.session_state.logged_in=True; st.rerun()
        else:st.error('Incorrect password')
    st.stop()

def api(method,path,**kwargs):
    if not BASE: raise RuntimeError('Configure WORKER_URL to enable QIIME 2 processing.')
    headers=kwargs.pop('headers',{}); headers['Authorization']='Bearer '+TOKEN
    r=requests.request(method,BASE+path,headers=headers,timeout=30,**kwargs)
    if not r.ok: raise RuntimeError(r.text[:1200])
    return r

def ppath(s=''):return '/p/'+st.session_state.pid+s

def params():return {'key':st.session_state.pkey}

def link(s):return PUBLIC+ppath(s)+'?'+urlencode(params())

def get(s):return api('GET',ppath(s),params=params()).json()

def choose(label,opts,key):
    result=st.selectbox(label,['— None —']+opts,key=key)
    return None if result=='— None —' else result

def download_table(d,label='Download CSV'):
    st.download_button(label,d.to_csv().encode(),'microbiome_table.csv','text/csv')

st.title('Microbiome Workbench')
st.caption('Oxford Nanopore full-length 16S · Emu + QIIME 2 · Data utilities')
page=st.sidebar.radio('Workspace',['Overview','1 · Project & upload','2 · Run analysis','3 · Results & logs','4 · Merge & stratify tables','5 · Visualisation'])
if st.sidebar.button('Sign out'):
    st.session_state.clear();st.rerun()
if st.session_state.get('pid'):st.sidebar.caption('Project '+st.session_state.pid[:10])
if page=='Overview':
    st.markdown('### From sequencing reads to microbiome results')
    for col,(num,title,desc) in zip(st.columns(4),[('01','Upload','FASTQ, QZA, metadata'),('02','Quality control','Read quality and primers'),('03','Analyse','Emu or DADA2, taxonomy, diversity'),('04','Explore','Tables, plots and exports')]):
        with col:st.metric(num,title);st.caption(desc)
    st.info('Use Project & upload to create a private project, then open the large-file uploader. The analysis server runs separately from Streamlit and continues a submitted job when the browser closes, provided the server stays running.')
    st.markdown('**Recommended workflow:** run quality control first → select platform-specific trimming → run denoising and taxonomy → choose rarefaction depth from the feature table summary → run downstream diversity using the generated artifacts.')
    st.warning('Only trusted research users should upload QIIME artifacts or classifiers. Do not use identifiable participant metadata. Demultiplexed, basecalled ONT full-length 16S FASTQ is supported. POD5 and shotgun metagenomes are not supported.')
    if not BASE:st.warning('Processing server is not configured. Table utilities and plots can still be used.')
    st.markdown('[QIIME 2 documentation](https://docs.qiime2.org/en/2026.10/) · [QIIME 2 View (local artifact viewing)](https://view.qiime2.org/)')
elif page=='1 · Project & upload':
    st.header('Project and inputs')
    try:
        if st.button('Create new project',type='primary'):
            result=api('POST','/projects').json();st.session_state.pid=result['id'];st.session_state.pkey=result['key'];st.rerun()
        with st.expander('Reopen an existing project'):
            pid=st.text_input('Project ID'); pkey=st.text_input('Project access key',type='password')
            if st.button('Open project'):
                api('GET','/p/'+pid+'/files',params={'key':pkey});st.session_state.pid=pid;st.session_state.pkey=pkey;st.rerun()
        if st.session_state.get('pid'):
            with st.expander('Save your project access details'):
                st.code(st.session_state.pid);st.code(st.session_state.pkey);st.caption('Keep these private. They allow access to all project inputs and results.')
            st.link_button('Open large-file uploader',link('/uploader'),type='primary');st.caption('100 GiB per file and 500 GiB input quota by default on a configured worker. Actual limits depend on worker settings and storage.')
            st.button('Refresh file list'); fs=get('/files')
            if fs:st.dataframe(pd.DataFrame([{'Filename':x['name'],'Size (MiB)':round(x['bytes']/1024**2,2)} for x in fs]),hide_index=True,use_container_width=True)
            else:st.info('No completed uploads yet.')
            with st.expander('Input formats and example manifest'):
                st.markdown('Upload demultiplexed Phred+33 FASTQ(.gz). For ONT, upload one file per sample or uniquely named chunks, with a manifest mapping all chunks to their sample ID. For other modes use one file per sample/read direction. Paired reads such as `Sample01_R1.fastq.gz` and `Sample01_R2.fastq.gz` are paired automatically. For unusual names, upload a tab-separated manifest with the exact uploaded filenames.')
                st.code('sample-id\tforward\treverse\nSample01\tSample01_R1.fastq.gz\tSample01_R2.fastq.gz');st.download_button('Paired manifest template','sample-id\tforward\treverse\nSample01\tSample01_R1.fastq.gz\tSample01_R2.fastq.gz\n','manifest.tsv')
                st.download_button('ONT multiple-chunk manifest template','sample-id\tfile\nSample01\tbarcode01_chunk001.fastq.gz\nSample01\tbarcode01_chunk002.fastq.gz\nSample02\tbarcode02_chunk001.fastq.gz\n','ont-manifest.tsv');st.download_button('Single/CCS manifest template','sample-id\tfile\nSample01\tSample01.fastq.gz\n','manifest-single.tsv');st.download_button('Metadata template','sample-id\tGroup\tVisit\nSample01\tControl\tBaseline\n','metadata.tsv')
            if fs:
                with st.expander('Remove an uploaded input'):
                    fname=st.selectbox('Filename',[x['name'] for x in fs]);confirm=st.checkbox('Confirm input removal')
                    if st.button('Delete input',disabled=not confirm):api('DELETE',ppath('/files/'+quote(fname,safe='')),params=params());st.rerun()
            with st.expander('Delete entire project'):
                ack=st.checkbox('Permanently delete all inputs, runs and results')
                if st.button('Delete project',disabled=not ack):api('DELETE',ppath(),params=params());del st.session_state.pid;del st.session_state.pkey;st.rerun()
    except Exception as e:st.error(str(e))
elif page=='2 · Run analysis':
    st.header('Configure and run QIIME 2')
    if not st.session_state.get('pid'):st.info('Create or reopen a project first.');st.stop()
    try: fs=get('/files')
    except Exception as e:st.error(str(e));st.stop()
    names=[f['name'] for f in fs]; qzas=[n for n in names if n.lower().endswith('.qza')]; txt=[n for n in names if n.lower().endswith(('.tsv','.csv','.txt'))]
    mode=st.radio('Starting material',['ont','paired','single','ccs','artifact'],format_func=lambda m:{'ont':'Oxford Nanopore full-length 16S FASTQ (recommended)','paired':'Illumina paired-end FASTQ','single':'Illumina single-end FASTQ','ccs':'PacBio HiFi / CCS full-length 16S','artifact':'Existing QIIME 2 feature tables / merge runs'}[m])
    cfg={'mode':mode,'stage':'full'}
    if mode=='ont':
        st.info('Basecalled, demultiplexed full-length 16S reads. For multiple FASTQ chunks per barcode/sample, select the manifest mapping sample IDs to unique uploaded filenames. Without it, each file is treated as a separate sample.')
        cfg['stage']=st.radio('Run type',['full','quality'],format_func=lambda s:'Automatic ONT profiling' if s=='full' else 'Quality report only')
        cfg['manifest']=choose('ONT manifest TSV (optional for one file per sample)',[n for n in txt if not n.endswith('.csv')],'ont-manifest')
        a,b,c=st.columns(3);cfg['min_len']=a.number_input('Minimum length (bp)',min_value=1,value=1200);cfg['max_len']=b.number_input('Maximum length (bp)',min_value=1,value=1800);cfg['min_q']=c.number_input('Minimum mean read Q score',min_value=0.0,max_value=60.0,value=10.0,step=1.0)
        st.caption('Starting thresholds for full-length 16S, not universal settings. Inspect kit, read lengths and control samples. Length filtering follows primer removal if selected.')
        if st.checkbox('Remove residual terminal primers (both required; discard reads without both)',value=False):
            a,b=st.columns(2);cfg['primer_f']=a.text_input('Forward primer (5′ → 3′)').strip();cfg['primer_r']=b.text_input('Reverse primer (5′ → 3′)').strip()
        st.caption('Optional Cutadapt handles both read orientations. Leave unchecked if sequencing provider already removed primers.')
        if cfg['stage']=='full':
            st.caption('Emu uses the reference database installed on the worker. Estimated counts can be fractional. No ASV inference, rarefaction or UniFrac in this mode.')
            cfg['levels']=[RANKS.index(x)+1 for x in st.multiselect('Export taxonomic levels',RANKS,default=RANKS)]
    elif mode!='artifact':
        cfg['stage']=st.radio('Run type',['quality','full'],format_func=lambda s:'Quality report only (recommended first)' if s=='quality' else 'Automatic pipeline with selected parameters')
        cfg['manifest']=choose('Custom manifest TSV (optional)',[n for n in txt if not n.endswith('.csv')],'manifest')
        if cfg['stage']=='full':
            st.subheader('Read processing');cfg['primer_f']=st.text_input('5′ forward primer (optional for Illumina; required for CCS)').strip();cfg['primer_r']=st.text_input('Reverse-read 5′ primer / CCS 3′ adapter in read orientation').strip() if mode!='single' else ''
            if mode=='ccs':
                a,b=st.columns(2);cfg['min_len']=a.number_input('Minimum read length',min_value=1,value=1000);cfg['max_len']=b.number_input('Maximum read length',min_value=1,value=1600);st.caption('Confirm the 3′ adapter orientation and expected amplicon length with your sequencing protocol.')
            else:
                a,b=st.columns(2);cfg['trunc_f']=a.number_input('Forward truncation length (0 = no truncation)',min_value=0,value=0);cfg['trim_f']=b.number_input('Forward trim-left bases',min_value=0,value=0)
                if mode=='paired':cfg['trunc_r']=a.number_input('Reverse truncation length (0 = no truncation)',min_value=0,value=0);cfg['trim_r']=b.number_input('Reverse trim-left bases',min_value=0,value=0);st.caption('Ensure sufficient forward/reverse overlap after trimming; never use universal truncation values.')
    else:
        cfg['tables']=st.multiselect('FeatureTable[Frequency] artifacts / BIOM 2.1 count tables',qzas+[n for n in names if n.lower().endswith('.biom')]);cfg['sequences']=st.multiselect('Representative sequence artifacts',qzas);cfg['taxonomies']=st.multiselect('Taxonomy artifacts',qzas);cfg['tree']=choose('Rooted tree artifact',qzas,'tree');cfg['sum_overlap']=st.checkbox('Sum counts for overlapping sample IDs (technical replicates only)')
        st.warning('Merge only compatible sequencing regions, feature IDs and preprocessing protocols. Use the same taxonomy reference. Do not sum biological replicates.')
    cfg['threads']=st.number_input('CPU threads',1,16,4)
    if mode=='ont':
        if cfg['stage']=='full':cfg['metadata']=choose('Sample metadata (optional)',txt,'ont-metadata')
        checked=st.checkbox('I confirm these are demultiplexed full-length 16S reads and have reviewed the parameters')
    elif cfg['stage']=='full':
        cfg['metadata']=choose('Sample metadata (required for artifact runs; optional for raw)',txt,'metadata')
        st.subheader('Taxonomy and diversity');cfg['classifier']=choose('Compatible trusted taxonomy classifier (.qza)',qzas,'classifier');st.caption('Upload a SILVA/other classifier trained for your sequenced region and QIIME/scikit-learn version. No classifier is bundled.')
        cfg['exclude_organelle']=st.checkbox('Remove mitochondrial and chloroplast reads')
        cfg['levels']=[RANKS.index(x)+1 for x in st.multiselect('Export taxonomic levels',RANKS,default=['Phylum','Family','Genus'])];cfg['phylogeny']=st.checkbox('Build phylogenetic tree (MAFFT + FastTree)',value=False);cfg['diversity']=st.checkbox('Run alpha/beta diversity and PCoA',value=False)
        if cfg['diversity']:
            cfg['depth']=st.number_input('Rarefaction depth (read counts/sample)',min_value=1,value=1000);cfg['group']=st.text_input('Optional categorical metadata column for PERMANOVA/PERMDISP');st.warning('Inspect sequencing depth and sample retention before choosing depth. Group tests assume independent samples; do not use ordinary PERMANOVA for repeated measurements without an appropriate design.')
        checked=st.checkbox('I have reviewed the quality/analysis parameters and input compatibility')
    else:checked=True
    if st.button('Start analysis',type='primary',disabled=not checked):
        try:
            if mode=='ont' and (cfg['min_len']>cfg['max_len'] or bool(cfg.get('primer_f'))!=bool(cfg.get('primer_r'))):raise ValueError('Check length bounds and supply both primers or neither.')
            if mode=='artifact' and (not cfg.get('tables') or not cfg.get('metadata')):raise ValueError('Select feature table(s) and metadata.')
            result=api('POST',ppath('/runs'),params=params(),json=cfg).json();st.success('Submitted run '+result['run'][:10]+'. Open Results & logs to track it.')
        except Exception as e:st.error(str(e))
elif page=='3 · Results & logs':
    st.header('Analysis results')
    if not st.session_state.get('pid'):st.info('Create or reopen a project first.');st.stop()
    try:
        st.button('Refresh status');runs=get('/runs')
        if not runs:st.info('No analysis has been submitted.');st.stop()
        rid=st.selectbox('Run',[r['id'] for r in runs],format_func=lambda x:next(r['state']+' · '+x[:10] for r in runs if r['id']==x));s=get('/runs/'+rid)
        st.metric('Status',s['state'].title());st.write(s['step']);st.code(s['log'] or 'No command output yet.',language='text')
        if s['state'] in ('running','queued'):st.info('Refresh to see current status. Do not download incomplete files.')
        else:
            st.markdown('### Output files');st.caption('Download QZV and open it in QIIME 2 View. Reuse QZA files directly in this project, or download them for another project. For confidential studies, use an approved local viewer.')
            if s['state']=='complete':
                eligible=[f['name'] for f in s['files'] if f['name'].endswith(('.qza','.tsv','.csv'))]
                if eligible:
                    with st.expander('Reuse a result in the next analysis (no re-upload)'):
                        source=st.selectbox('Result to reuse',eligible); dest=st.text_input('New input filename',value=rid[:8]+'-'+source.split('/')[-1],key='reuse-'+rid+'-'+source)
                        if st.button('Add to project inputs'):
                            api('POST',ppath('/runs/'+rid+'/reuse'),params=params(),json={'source':source,'name':dest});st.success('Result added to project inputs.')
            for f in s['files']:
                a,b,c=st.columns([5,1,2]);a.write(f['name']);b.caption(f"{f['bytes']/1024**2:.2f} MiB");c.link_button('Download',link('/runs/'+rid+'/download/'+quote(f['name'],safe='/')))
    except Exception as e:st.error(str(e))
elif page=='4 · Merge & stratify tables':
    st.header('Merge and stratify abundance tables');st.caption('CSV/TSV: first column feature or taxon ID; remaining columns sample counts or abundances. Never sum overlapping relative-abundance columns; use compatible estimated counts for technical-replicate merges only. Use QIIME artifact mode for large BIOM/QZA datasets.')
    uploads=st.file_uploader('Count tables',type=['csv','tsv','txt'],accept_multiple_files=True); operation=st.selectbox('Operation',['Merge tables by feature IDs','Collapse by taxonomic rank','Collapse each table then merge','Stratify / aggregate by metadata group'])
    overlap=st.radio('Duplicate sample IDs',['error','sum'],format_func=lambda x:'Reject duplicates (recommended)' if x=='error' else 'Sum technical replicates only',horizontal=True)
    taxfile=st.file_uploader('Taxonomy TSV (Feature ID, Taxon)',type=['tsv','txt','csv']) if 'ollapse' in operation else None
    rank=st.selectbox('Taxonomic rank',RANKS,index=5) if taxfile else None
    meta=st.file_uploader('Metadata CSV / TSV (first column sample ID)',type=['csv','tsv']) if 'metadata' in operation else None
    try:
        m=None;group=None;method=None
        if meta:
            m=pd.read_csv(meta,sep='\t' if meta.name.endswith('.tsv') else ',',index_col=0,dtype=str);m.index=m.index.astype(str);group=st.selectbox('Stratify by',m.columns);method=st.selectbox('Aggregation',['sum','mean','median'])
        if st.button('Process tables',type='primary'):
            ds=[load_counts(f) for f in uploads]
            if not ds:raise ValueError('Upload at least one table.')
            if 'ollapse' in operation:
                if not taxfile:raise ValueError('Upload taxonomy.')
                taxdf=pd.read_csv(taxfile,sep=',' if taxfile.name.endswith('.csv') else '\t',index_col=0,dtype=str);taxdf.index=taxdf.index.astype(str)
                if taxdf.index.duplicated().any() or 'Taxon' not in taxdf:raise ValueError('Taxonomy must have unique feature IDs and a Taxon column.')
                level=RANKS.index(rank)+1
                d=merge_counts([collapse(t,taxdf['Taxon'],level) for t in ds],overlap) if operation=='Collapse each table then merge' else collapse(merge_counts(ds,overlap),taxdf['Taxon'],level)
            else:
                d=merge_counts(ds,overlap)
                if 'metadata' in operation:
                    if m is None:raise ValueError('Upload metadata.')
                    d=aggregate(d,m,group,method)
            st.session_state.table_result=d
        if 'table_result' in st.session_state:
            d=st.session_state.table_result;st.success(f'{d.shape[0]:,} features × {d.shape[1]:,} samples/groups');st.dataframe(d.head(100));download_table(d)
    except Exception as e:st.error(str(e))
elif page=='5 · Visualisation':
    st.header('Explore microbiome abundance');f=st.file_uploader('Count table CSV / TSV',type=['csv','tsv','txt']);use=st.checkbox('Use latest table utility result',value=False,disabled='table_result' not in st.session_state)
    try:
        d=st.session_state.table_result if use else load_counts(f) if f else None
        if d is None:st.info('Upload a count table or process one in Merge & stratify tables.');st.stop()
        plot=st.selectbox('Plot',['Stacked relative abundance','Abundance heatmap','Richness / Shannon boxplot','Table column totals','Top taxa bar chart','Bray–Curtis PCoA'])
        top=st.slider('Top taxa',5,100,20);n=min(top,len(d));r=relative(d);keep=r.mean(axis=1).nlargest(n).index
        if plot=='Stacked relative abundance':
            x=r.loc[keep].copy();other=r.drop(index=keep).sum(axis=0,min_count=1)
            if len(keep)<len(r):x.loc['Other (remaining taxa)']=other
            fig=px.bar(x.T,barmode='stack',labels={'index':'Sample','value':'Relative abundance (%)','variable':'Taxon'})
        elif plot=='Abundance heatmap':fig=px.imshow(r.loc[keep],aspect='auto',labels={'color':'Relative abundance (%)'})
        elif plot=='Table column totals':fig=px.bar(d.sum().rename('Total'),labels={'index':'Sample','value':'Input table total'});st.caption('Not sequencing depth for relative-abundance or filtered estimated-count tables; use ONT QC for observed reads.')
        elif plot=='Top taxa bar chart':fig=px.bar(r.mean(axis=1).nlargest(n).sort_values().rename('Mean abundance (%)'),orientation='h',labels={'index':'Taxon','value':'Mean abundance (%)'})
        elif plot=='Richness / Shannon boxplot':
            p=r/100;metrics=pd.DataFrame({'Observed features':(d>0).sum(),'Shannon':-(p.where(p>0)*np.log(p.where(p>0))).sum(min_count=1)});metric=st.selectbox('Metric',metrics.columns);fig=px.box(metrics,y=metric,points='all');st.warning('Descriptive, unrarefied metrics. ONT taxon profiles measure observed taxa, not ASVs; do not apply count-based rarefaction to Emu estimates.')
        else:
            from scipy.spatial.distance import pdist,squareform
            if d.shape[1]>1500:raise ValueError('For >1,500 samples use the QIIME 2 worker for ordination.')
            x=d.loc[:,d.sum()>0].T
            if len(x)<3:raise ValueError('At least three nonempty samples required.')
            D=squareform(pdist(x.div(x.sum(axis=1),axis=0),metric='braycurtis'));J=np.eye(len(x))-np.ones((len(x),len(x)))/len(x);vals,vecs=np.linalg.eigh(-.5*J@(D**2)@J);order=np.argsort(vals)[::-1];vals=vals[order];vecs=vecs[:,order];pos=vals[vals>0].sum();coords=vecs[:,:2]*np.sqrt(np.maximum(vals[:2],0));df=pd.DataFrame(coords,index=x.index,columns=['PCo1','PCo2']);df['Sample']=df.index
            fig=px.scatter(df,x='PCo1',y='PCo2',hover_name='Sample',labels={'PCo1':f'PCo1 ({100*max(vals[0],0)/pos:.1f}%)' if pos else 'PCo1','PCo2':f'PCo2 ({100*max(vals[1],0)/pos:.1f}%)' if pos else 'PCo2'});st.caption('Bray–Curtis on relative abundances; percentages use positive eigenvalues. Negative eigenvalues may indicate non-Euclidean distances.')
        fig.update_layout(template='plotly_white',height=580);st.plotly_chart(fig,use_container_width=True);st.download_button('Download interactive HTML',fig.to_html(include_plotlyjs=True),'microbiome_plot.html','text/html')
    except Exception as e:st.error(str(e))
