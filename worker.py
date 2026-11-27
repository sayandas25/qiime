"""Private analysis worker; one queued job at a time, persistent disk."""
import os, re, json, secrets, shutil, threading, hmac, uuid
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from fastapi import FastAPI, HTTPException, Request, Header, Query
from fastapi.responses import HTMLResponse, FileResponse, PlainTextResponse
from pipeline import execute

ROOT=Path(os.environ.get('PROJECT_ROOT','/data/projects')).resolve(); ROOT.mkdir(parents=True,exist_ok=True)
if os.environ.get('TMPDIR'): Path(os.environ['TMPDIR']).mkdir(parents=True,exist_ok=True)
TOKEN=os.environ.get('WORKER_TOKEN',''); LIMIT=int(os.environ.get('MAX_FILE_GIB','100'))*1024**3; QUOTA=int(os.environ.get('PROJECT_GIB','500'))*1024**3
CHUNK=16*1024**2; pool=ThreadPoolExecutor(max_workers=1); lock=threading.RLock(); app=FastAPI(docs_url=None,redoc_url=None,openapi_url=None)

def fail(msg,code=400): raise HTTPException(code,msg)
def write(p,v):
    tmp=p.with_suffix('.tmp'); tmp.write_text(json.dumps(v,indent=2)); tmp.replace(p)
def read(p): return json.loads(p.read_text())
def safe(name):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,199}',name) or name in ('.','..'): fail('Use safe filenames (letters, numbers, dots, underscores, hyphens).')
    return name
def project(pid,key):
    if not re.fullmatch('[a-f0-9]{32}',pid): fail('Unknown project',404)
    p=ROOT/pid
    if not (p/'access.json').is_file() or not hmac.compare_digest(read(p/'access.json')['key'],key or ''): fail('Invalid project access',403)
    return p
def runpath(p,rid):
    if not re.fullmatch('[a-f0-9]{32}',rid) or not (p/'runs'/rid/'status.json').exists(): fail('Unknown run',404)
    return p/'runs'/rid
def busy(p): return any(read(f)['state'] in ('queued','running') for f in (p/'runs').glob('*/status.json'))
# A stopped worker cannot resume an interrupted external command automatically.
for status in ROOT.glob('*/runs/*/status.json'):
    s=read(status)
    if s['state'] in ('queued','running'): write(status,dict(s,state='failed',step='Worker restarted; start a new run.'))

@app.post('/projects')
def create(authorization:str=Header(default='')):
    if not TOKEN or not hmac.compare_digest(authorization,'Bearer '+TOKEN): fail('Worker authorization required',403)
    pid=uuid.uuid4().hex; key=secrets.token_urlsafe(32); p=ROOT/pid; p.mkdir(); (p/'inputs').mkdir(); (p/'partial').mkdir(); (p/'runs').mkdir(); write(p/'access.json',{'key':key}); return {'id':pid,'key':key}

@app.get('/p/{pid}/files')
def files(pid:str,key:str):
    p=project(pid,key); return [{'name':f.name,'bytes':f.stat().st_size} for f in sorted((p/'inputs').iterdir()) if f.is_file()]

@app.get('/p/{pid}/offset/{name}')
def offset(pid:str,name:str,key:str):
    p=project(pid,key); safe(name); partial=p/'partial'/name; return {'offset':partial.stat().st_size if partial.exists() else 0,'complete':(p/'inputs'/name).exists()}

@app.put('/p/{pid}/upload/{name}')
async def upload(pid:str,name:str,request:Request,key:str,offset:int=Query(ge=0),total:int=Query(gt=0)):
    p=project(pid,key); safe(name)
    if total>LIMIT: fail('File exceeds configured per-file limit',413)
    length=int(request.headers.get('content-length','0'))
    if not 0<length<=CHUNK or offset+length>total: fail('Invalid chunk length',413)
    buf=bytearray()
    async for block in request.stream():
        buf.extend(block)
        if len(buf)>length: fail('Chunk exceeds declared length',413)
    data=buf
    if len(data)!=length: fail('Incomplete chunk')
    with lock:
        if busy(p): fail('Cannot modify inputs during analysis',409)
        dst=p/'inputs'/name; part=p/'partial'/name
        if dst.exists(): fail('File already uploaded; use a different filename',409)
        current=part.stat().st_size if part.exists() else 0
        if offset!=current: fail('Offset mismatch; reselect file to resume',409)
        # Persist total for resume; selecting a different same-name/same-size file is not safe.
        info=p/'partial'/('.'+name+'.json')
        if info.exists() and read(info)['total']!=total: fail('Different file size; reset partial upload first',409)
        usage=sum(f.stat().st_size for d in ('inputs','partial') for f in (p/d).iterdir() if f.is_file())
        if usage+length>QUOTA or shutil.disk_usage(ROOT).free<length+1024**3: fail('Storage quota or free space exceeded',507)
        write(info,{'total':total});
        with part.open('ab') as f: f.write(data)
        if offset+length==total: part.replace(dst); info.unlink(missing_ok=True)
    return {'offset':offset+length,'complete':offset+length==total}

@app.delete('/p/{pid}/partial/{name}')
def reset(pid:str,name:str,key:str):
    p=project(pid,key); safe(name)
    with lock:
        if busy(p): fail('Analysis active',409)
        (p/'partial'/name).unlink(missing_ok=True); (p/'partial'/('.'+name+'.json')).unlink(missing_ok=True)
    return {'ok':True}

@app.delete('/p/{pid}/files/{name}')
def delete_input(pid:str,name:str,key:str):
    p=project(pid,key); safe(name)
    with lock:
        if busy(p): fail('Analysis active',409)
        (p/'inputs'/name).unlink(missing_ok=True)
    return {'ok':True}

@app.get('/p/{pid}/uploader',response_class=HTMLResponse)
def uploader(pid:str,key:str):
    project(pid,key); return HTMLResponse(Path(__file__).with_name('uploader.html').read_text(),headers={'Referrer-Policy':'no-referrer','Cache-Control':'no-store'})

def validate(cfg,p):
    if cfg.get('mode') not in ('ont','paired','single','ccs','artifact') or cfg.get('stage') not in ('quality','full'): fail('Invalid workflow')
    if cfg['mode']=='artifact' and cfg['stage']=='quality': fail('QC is for raw data')
    for field in ('metadata','manifest','classifier','tree'):
        if cfg.get(field): safe(cfg[field]);
        if cfg.get(field) and not (p/'inputs'/cfg[field]).is_file(): fail('Missing '+field)
    for field in ('tables','sequences','taxonomies'):
        for name in cfg.get(field,[]):
            safe(name)
            if not (p/'inputs'/name).is_file(): fail('Missing '+name)
    for field in ('threads','trunc_f','trunc_r','trim_f','trim_r','min_len','max_len','depth'):
        if field in cfg and (type(cfg[field]) is not int or not 0<=cfg[field]<=1000000000): fail('Invalid '+field)
    if not 1<=cfg.get('threads',4)<=int(os.environ.get('MAX_THREADS','16')): fail('Threads exceed worker limit')
    if any(not isinstance(x,int) or not 1<=x<=7 for x in cfg.get('levels',[])): fail('Invalid taxonomic rank')
    for field in ('primer_f','primer_r'):
        if cfg.get(field) and not re.fullmatch('[ACGTRYSWKMBDHVNacgtryswkmbdhvn]{1,200}',cfg[field]): fail('Invalid primer sequence')
    if cfg['mode']=='ont':
        if not 1<=cfg.get('min_len',1200)<=cfg.get('max_len',1800): fail('Invalid ONT length bounds')
        if type(cfg.get('min_q',10)) not in (int,float) or not 0<=cfg.get('min_q',10)<=60: fail('Invalid ONT quality threshold')
        if bool(cfg.get('primer_f')) != bool(cfg.get('primer_r')): fail('Supply both primers or neither')
    if cfg.get('diversity') and cfg.get('depth',0)<1: fail('Set rarefaction depth')

def process(p,r,cfg):
    s=r/'status.json'
    def report(step): write(s,{'state':'running','step':step})
    try:
        report('Starting QIIME 2'); execute(p,r,cfg,report); write(s,{'state':'complete','step':'Finished'})
    except Exception as e:
        with (r/'commands.log').open('a') as f: f.write('\nERROR: '+str(e)+'\n')
        write(s,{'state':'failed','step':str(e)})

@app.post('/p/{pid}/runs')
def start(pid:str,cfg:dict,key:str):
    p=project(pid,key); validate(cfg,p)
    with lock:
        if busy(p): fail('This project already has an active analysis',409)
        rid=uuid.uuid4().hex; r=p/'runs'/rid; r.mkdir(); write(r/'parameters.json',cfg); write(r/'status.json',{'state':'queued','step':'Waiting for worker'}); pool.submit(process,p,r,cfg)
    return {'run':rid}

@app.get('/p/{pid}/runs')
def runs(pid:str,key:str):
    p=project(pid,key); return [dict(id=r.name,**read(r/'status.json')) for r in sorted((p/'runs').iterdir(),key=lambda x:x.stat().st_mtime,reverse=True)]

@app.get('/p/{pid}/runs/{rid}')
def status(pid:str,rid:str,key:str):
    r=runpath(project(pid,key),rid); log=r/'commands.log'; tail=''
    if log.exists():
        with log.open('rb') as f: f.seek(max(0,log.stat().st_size-20000)); tail=f.read().decode(errors='replace')
    return dict(**read(r/'status.json'),log=tail,files=[{'name':str(f.relative_to(r)),'bytes':f.stat().st_size} for f in sorted(r.rglob('*')) if f.is_file() and not f.name.endswith('.tmp')])

@app.post('/p/{pid}/runs/{rid}/reuse')
def reuse(pid:str,rid:str,cfg:dict,key:str):
    p=project(pid,key); r=runpath(p,rid); source=cfg.get('source',''); name=safe(cfg.get('name',''))
    f=(r/source).resolve()
    if not f.is_relative_to(r) or not f.is_file() or f.suffix not in ('.qza','.tsv','.csv'): fail('Select an artifact or table output')
    if read(r/'status.json')['state']!='complete': fail('Run must be complete')
    with lock:
        if busy(p): fail('Analysis active',409)
        dest=p/'inputs'/name
        if dest.exists() or (p/'partial'/name).exists(): fail('Input filename already in use',409)
        if sum(x.stat().st_size for d in ('inputs','partial') for x in (p/d).iterdir() if x.is_file())+f.stat().st_size>QUOTA: fail('Input quota exceeded',507)
        os.link(f,dest)
    return {'ok':True}

@app.get('/p/{pid}/runs/{rid}/download/{name:path}')
def download(pid:str,rid:str,name:str,key:str):
    r=runpath(project(pid,key),rid); f=(r/name).resolve()
    if not f.is_relative_to(r) or not f.is_file(): fail('File not found',404)
    return FileResponse(f,filename=f.name,headers={'Cache-Control':'no-store','Referrer-Policy':'no-referrer'})

@app.delete('/p/{pid}')
def remove(pid:str,key:str):
    p=project(pid,key)
    with lock:
        if busy(p): fail('Analysis active',409)
        shutil.rmtree(p)
    return {'ok':True}
