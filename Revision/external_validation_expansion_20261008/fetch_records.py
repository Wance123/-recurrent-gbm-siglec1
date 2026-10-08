
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import urllib.request,json,hashlib
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).parent
def fetch(acc):
 url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full'
 try:
  with urllib.request.urlopen(url,timeout=45) as r:b=r.read(500001)
  if len(b)>500000 or not b.startswith(b'^'):raise ValueError('not bounded GEO text record')
  (ROOT/'input'/(acc+'_record.txt')).write_bytes(b)
  return dict(accession=acc,url=url,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),status='downloaded')
 except Exception as e:return dict(accession=acc,url=url,status='error',error=str(e))
if __name__=='__main__':
 with ThreadPoolExecutor(max_workers=3) as p:r=list(p.map(fetch,['GSE212067','GSE42670','GSE7696']))
 (ROOT/'record_download_audit.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
 print(json.dumps(r,indent=2),flush=True)
