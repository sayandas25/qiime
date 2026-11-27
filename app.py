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
PASSWORD=setting('PORTAL_PASSWORD')
if not PASSWORD:
    st.error('Set PORTAL_PASSWORD in Streamlit secrets or environment before using this private portal.'); st.stop()
if not st.session_state.get('logged_in'):
    st.title('Microbiome Workbench'); st.caption('Private research workspace · Galaxy Europe')
    with st.form('login'):
        pwd=st.text_input('Workspace password',type='password'); submit=st.form_submit_button('Sign in')
    if submit:
        if hmac.compare_digest(pwd,PASSWORD): st.session_state.logged_in=True; st.rerun()
        else:st.error('Incorrect password')
    st.stop()

from bioblend.galaxy import GalaxyInstance
GALAXY='https://usegalaxy.eu'
TRAINING='https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/nanopore-16S-metagenomics/workflows/main_workflow.ga'

def download_table(d,label='Download CSV'):
    st.download_button(label,d.to_csv().encode(),'microbiome_table.csv','text/csv')

def galaxy():
    key=st.session_state.get('galaxy_key','')
    if not key:
        st.info('Enter your personal Galaxy API key in the sidebar.');st.stop()
    return GalaxyInstance(url=GALAXY,key=key)

def safe_error():
    st.error('Galaxy request failed. Check the API key, workflow inputs, account quota and Galaxy status. See Galaxy for detailed errors.')

def pick(label,items,key):
    ids=[x['id'] for x in items]
    labels={x['id']:f"{x.get('name','Unnamed')} · {x.get('state','')} · {x['id'][:8]}" for x in items}
    return st.selectbox(label,ids,index=None,format_func=lambda i:labels[i],key=key)

def history(gi,key):
    hs=gi.histories.get_histories(deleted=False)
    return pick('Galaxy history',hs,key)

def content(gi,hid):
    return [x for x in gi.histories.show_history(hid,contents=True,deleted=False) if not x.get('deleted') and not x.get('purged')]

st.title('Microbiome Workbench')
st.caption('Oxford Nanopore full-length 16S · Galaxy Europe · Data utilities')
page=st.sidebar.radio('Workspace',['Overview','1 · Project & upload','2 · Run analysis','3 · Results & logs','4 · Merge & stratify tables','5 · Visualisation'])
if st.sidebar.button('Sign out'):
    st.session_state.clear();st.rerun()
st.sidebar.text_input('Your Galaxy API key',type='password',key='galaxy_key',help='User → Preferences → Manage API key. Stored in this browser session on the Streamlit server; do not share or commit.')
st.sidebar.link_button('Open Galaxy Europe',GALAXY)
if page=='Overview':
    st.markdown('### From Nanopore reads to microbiome results')
    st.info('Upload large FASTQ files directly into your Galaxy history. The portal submits an imported workflow to Galaxy and monitors its outputs. No private worker server is needed.')
    st.markdown('1. Connect your Galaxy API key.\n2. Create a history and upload files directly in Galaxy.\n3. Import and review the Nanopore 16S workflow; create a list collection if required.\n4. Select the history and workflow here; submit and track results.\n5. Download compatible abundance tables and explore them in the data utilities.')
    st.link_button('Official Nanopore 16S tutorial','https://training.galaxyproject.org/training-material/topics/microbiome/tutorials/nanopore-16S-metagenomics/tutorial.html')
    st.code(TRAINING,language=None)
    st.warning('The training workflow is a Kraken2/SILVA example, not an Emu or QIIME 2 ASV pipeline. Review primers, QC, database version and controls before use. Use public test data first; do not upload sensitive human-derived data without approvals and compliance with Galaxy terms. Galaxy storage quotas and queue limits apply.')
elif page=='1 · Project & upload':
    st.header('Project and inputs')
    gi=galaxy()
    try:
        name=st.text_input('New history name','Nanopore 16S project')
        if st.button('Create Galaxy history',type='primary'):
            if not name.strip():raise ValueError('Enter a history name')
            gi.histories.create_history(name=name.strip());st.success('History created.')
        st.link_button('Upload FASTQ directly to Galaxy',GALAXY,type='primary')
        st.caption('In Galaxy select your history, click Upload, choose files and start. For workflow collection inputs, build a list collection with one sample per element; concatenate chunks belonging to the same sample in Galaxy first. No raw FASTQ passes through Streamlit.')
        st.button('Refresh histories and files')
        hid=history(gi,'upload_history')
        if hid:
            rows=content(gi,hid)
            if rows:st.dataframe(pd.DataFrame([{k:x.get(k) for k in ['name','state','history_content_type','hid','id']} for x in rows]),hide_index=True)
            else:st.info('No datasets yet.')
    except Exception:safe_error()
elif page=='2 · Run analysis':
    st.header('Submit a Galaxy workflow')
    st.link_button('Import / review workflow in Galaxy',GALAXY+'/workflows/list')
    st.caption('Import this official training workflow URL in Galaxy, then edit and validate settings for your own 16S sequencing kit. This portal runs the saved settings; it does not automatically tune them.')
    st.code(TRAINING,language=None)
    gi=galaxy()
    try:
        st.button('Refresh workflows')
        hid=history(gi,'run_history')
        wfs=gi.workflows.get_workflows()
        wid=pick('Imported workflow',wfs,'workflow_id')
        if not hid or not wid:st.stop()
        wf=gi.workflows.show_workflow(wid)
        inputs=wf.get('inputs',{})
        if not inputs:st.warning('No dataset inputs reported by Galaxy. Run this workflow in Galaxy directly.');st.stop()
        rows=content(gi,hid)
        singles=[x for x in rows if x.get('history_content_type')=='dataset' and x.get('state')=='ok']
        collections=[x for x in rows if x.get('history_content_type')=='dataset_collection']
        selected={};valid=True
        for idx,info in inputs.items():
            step=wf.get('steps',{}).get(str(idx),{})
            typ=step.get('type','') or info.get('type','')
            if typ in ('parameter_input','parameter'):
                st.warning('This workflow has a runtime parameter input; configure and run it in Galaxy directly.');valid=False;continue
            if typ in ('data_collection_input','dataset_collection'):
                options=collections;kind='hdca'
            elif typ in ('data_input','dataset'):
                options=singles;kind='hda'
            else:
                kind=st.radio('Input type: '+str(info.get('label',idx)),['hda','hdca'],format_func=lambda x:'Dataset' if x=='hda' else 'Collection',key='kind'+wid+str(idx))
                options=singles if kind=='hda' else collections
            val=pick(str(info.get('label',idx)),options,'input'+hid+wid+str(idx)+kind)
            if val:selected[str(idx)]={'id':val,'src':kind}
            else:valid=False
        ack=st.checkbox('I reviewed the saved workflow settings, database and input sample mapping; data are permitted on Galaxy.')
        if st.button('Start workflow',type='primary',disabled=not (valid and ack)):
            result=gi.workflows.invoke_workflow(wid,inputs=selected,history_id=hid,inputs_by='step_index',require_exact_tool_versions=True)
            st.success('Submitted to Galaxy. Track the jobs and output datasets in Results & logs.')
            st.code(str(result.get('id','See Galaxy history')),language=None)
    except Exception:safe_error()
elif page=='3 · Results & logs':
    st.header('Galaxy results and job status')
    gi=galaxy()
    try:
        st.button('Refresh job status')
        hid=history(gi,'result_history')
        if not hid:st.stop()
        st.link_button('Open full history in Galaxy',GALAXY+'/histories/view?id='+hid)
        rows=content(gi,hid)
        if rows:st.dataframe(pd.DataFrame([{k:x.get(k) for k in ['name','state','history_content_type','hid','id']} for x in rows]),hide_index=True)
        jobs=gi.jobs.get_jobs(history_id=hid)
        if jobs:st.dataframe(pd.DataFrame([{k:x.get(k) for k in ['tool_id','state','create_time','update_time']} for x in jobs]),hide_index=True)
        st.caption('Download outputs from Galaxy in your signed-in browser, including collection elements, FASTQ, reports and QIIME artifacts when produced. The portal does not proxy large results or automatically convert Kraken reports to abundance matrices.')
    except Exception:safe_error()
elif page=='4 · Merge & stratify tables':
    st.header('Merge and stratify abundance tables');st.caption('CSV/TSV: first column feature or taxon ID; remaining columns sample counts or abundances. Never sum overlapping relative-abundance columns; use compatible estimated counts for technical-replicate merges only. Use Galaxy for large BIOM/QZA datasets.')
    uploads=st.file_uploader('Abundance tables',type=['csv','tsv','txt'],accept_multiple_files=True); operation=st.selectbox('Operation',['Merge tables by feature IDs','Collapse by taxonomic rank','Collapse each table then merge','Stratify / aggregate by metadata group'])
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
            if d.shape[1]>1500:raise ValueError('For >1,500 samples use a Galaxy workflow for ordination.')
            x=d.loc[:,d.sum()>0].T
            if len(x)<3:raise ValueError('At least three nonempty samples required.')
            D=squareform(pdist(x.div(x.sum(axis=1),axis=0),metric='braycurtis'));J=np.eye(len(x))-np.ones((len(x),len(x)))/len(x);vals,vecs=np.linalg.eigh(-.5*J@(D**2)@J);order=np.argsort(vals)[::-1];vals=vals[order];vecs=vecs[:,order];pos=vals[vals>0].sum();coords=vecs[:,:2]*np.sqrt(np.maximum(vals[:2],0));df=pd.DataFrame(coords,index=x.index,columns=['PCo1','PCo2']);df['Sample']=df.index
            fig=px.scatter(df,x='PCo1',y='PCo2',hover_name='Sample',labels={'PCo1':f'PCo1 ({100*max(vals[0],0)/pos:.1f}%)' if pos else 'PCo1','PCo2':f'PCo2 ({100*max(vals[1],0)/pos:.1f}%)' if pos else 'PCo2'});st.caption('Bray–Curtis on relative abundances; percentages use positive eigenvalues. Negative eigenvalues may indicate non-Euclidean distances.')
        fig.update_layout(template='plotly_white',height=580);st.plotly_chart(fig,use_container_width=True);st.download_button('Download interactive HTML',fig.to_html(include_plotlyjs=True),'microbiome_plot.html','text/html')
    except Exception as e:st.error(str(e))
