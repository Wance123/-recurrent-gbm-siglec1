
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import gzip, hashlib, json, urllib.request, re
from concurrent.futures import ThreadPoolExecutor

ROOT=Path(__file__).parent
OUT=ROOT/'input'; OUT.mkdir(exist_ok=True)
SERIES=['GSE274546','GSE212067','GSE4271','GSE42670','GSE7696','GSE79671']

def fetch(acc):
    bucket=acc[:-3]+'nnn'
    url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{acc}/soft/{acc}_family.soft.gz'
    p=OUT/(acc+'_family.soft.gz')
    try:
        if not p.exists():
            with urllib.request.urlopen(url,timeout=45) as r:
                b=r.read(45_000_001)
            if len(b)>45_000_000: raise ValueError('metadata download exceeds cap')
            gzip.decompress(b)
            p.write_bytes(b)
        b=p.read_bytes(); text=gzip.decompress(b).decode('utf-8',errors='replace')
        samples=[]; section=''; row=None; series=[]; in_table=False
        for line in text.splitlines():
            if line.startswith('^'):
                if row: samples.append(row)
                row=None; section=line.split(' = ')[0][1:]
                if section=='SAMPLE': row={'accession':line.split(' = ')[1]}
            if row is not None and line.startswith('!Sample_') and ' = ' in line:
                k,v=line[8:].split(' = ',1);row.setdefault(k,[]).append(v)
            if line.startswith('!Series_'): series.append(line)
        if row:samples.append(row)
        result={'accession':acc,'source_url':url,'record_url':f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}',
                'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'series':series,'samples':samples}
        (OUT/(acc+'_metadata.json')).write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
        print(acc,'samples',len(samples),'bytes',len(b),flush=True)
        return {'accession':acc,'status':'downloaded','samples':len(samples),'bytes':len(b)}
    except Exception as e:
        print(acc,type(e).__name__,str(e),flush=True)
        return {'accession':acc,'status':'error','error':str(e)}

if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=3) as pool: results=list(pool.map(fetch,SERIES))
    (ROOT/'metadata_download_audit.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
