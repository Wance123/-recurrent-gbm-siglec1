
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import urllib.request,re,gzip,json,hashlib
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).parent
def download(url,name,cap=45_000_000):
 p=ROOT/'input'/name
 if not p.exists():
  with urllib.request.urlopen(url,timeout=60) as r:b=r.read(cap+1)
  if len(b)>cap:raise ValueError('processed file exceeds download cap')
  gzip.decompress(b);p.write_bytes(b)
 b=p.read_bytes();return dict(file=name,url=url,bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
def matrix(acc):
 base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/matrix/'
 with urllib.request.urlopen(base,timeout=45) as r:t=r.read(100000).decode()
 names=sorted(set(re.findall(r'href="([^"]+series_matrix[^"/]*\.gz)"',t)))
 if not names:raise ValueError('no matrix file')
 return [download(base+n,n) for n in names]
def platform(acc):
 try:
  url=f'https://ftp.ncbi.nlm.nih.gov/geo/platforms/{acc[:-3]}nnn/{acc}/annot/{acc}.annot.gz'
  return [download(url,acc+'.annot.gz')]
 except Exception as e:print(acc,'annotation route',str(e),flush=True)
 url=f'https://ftp.ncbi.nlm.nih.gov/geo/platforms/{acc[:-3]}nnn/{acc}/soft/{acc}_family.soft.gz'
 return [download(url,acc+'_family.soft.gz',150_000_000)]
def run(kind,acc):
 try:
  r=(matrix if kind=='matrix' else platform)(acc);print(acc,[x['bytes'] for x in r],flush=True);return dict(accession=acc,status='downloaded',files=r)
 except Exception as e:print(acc,str(e),flush=True);return dict(accession=acc,status='error',error=str(e))
if __name__=='__main__':
 jobs=[('matrix','GSE212067'),('matrix','GSE7696'),('platform','GPL23126'),('platform','GPL570')]
 with ThreadPoolExecutor(max_workers=3) as p:r=list(p.map(lambda a:run(*a),jobs))
 (ROOT/'matrix_download_audit.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
