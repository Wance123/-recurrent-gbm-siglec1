
# Public release: path relocation only; analytical operations are unchanged.
import sys as _sys
from pathlib import Path as _ReleasePath
_release_root=next(p for p in _ReleasePath(__file__).resolve().parents if (p/'portable_paths.py').is_file())
_sys.path.insert(0,str(_release_root))
from portable_paths import project_path
from pathlib import Path
import csv,gzip,json,re,urllib.request,hashlib,io
import numpy as np
import pandas as pd
import evaluate_external as ev
ROOT=Path(__file__).parent;IN=ROOT/'input'
GENES=ev.GENES

def mapping_row(fields,header):
    h={k.lower():i for i,k in enumerate(header)}
    if 'gene symbol' in h:
        symbols=set(fields[h['gene symbol']].split(' /// '))
    elif 'gene_assignment' in h:
        symbols=set()
        for item in fields[h['gene_assignment']].split(' /// '):
            bits=item.split(' // ')
            if len(bits)>=2:symbols.add(bits[1].strip())
    else:raise ValueError('unrecognized official gene mapping field')
    symbols.discard('');symbols.discard('---')
    return next(iter(symbols)) if len(symbols)==1 and symbols<=set(GENES) else None

def fetch_clariom_mapping():
    dest=IN/'GPL23126_target_probe_mapping.json'
    if dest.exists():return json.loads(dest.read_text())['mapping']
    url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL23126&targ=self&form=text&view=full'
    mapping={};sha=hashlib.sha256();count=0;inside=False;complete=False;header=None
    cache=IN/'GPL23126_self_record.gz'
    if not cache.exists():
        request=urllib.request.Request(url,headers={'Accept-Encoding':'gzip'})
        with urllib.request.urlopen(request,timeout=60) as response:
            compressed=response.headers.get('Content-Encoding')=='gzip'
            b=response.read(150_000_001)
        if len(b)>150_000_000:raise ValueError('platform compressed transfer exceeds cap')
        cache.write_bytes(b if compressed else gzip.compress(b))
    # Stream the already bounded, complete official annotation; avoid materializing
    # its very large decompressed GO/transcript fields or downloading it again.
    with gzip.open(cache,'rb') as f:
        for b in f:
            count+=len(b);sha.update(b)
            line=b.decode('utf-8',errors='replace').strip('\r\n')
            if line.startswith('!platform_table_begin'):inside=True;continue
            if line.startswith('!platform_table_end'):complete=True;break
            if not inside or not line:continue
            v=next(csv.reader([line],delimiter='\t'))
            if header is None:header=v;continue
            gene=mapping_row(v,header)
            if gene:mapping[v[0]]=gene
    assert complete and 'SIGLEC1' in mapping.values()
    dest.write_text(json.dumps(dict(source_url=url,bytes_read=count,sha256_stream_to_table_end=sha.hexdigest(),complete_platform_table=True,mapping=mapping),indent=2),encoding='utf-8')
    print('Clariom probes',pd.Series(mapping).value_counts().to_dict(),flush=True)
    return mapping

def annot_mapping(path):
    mapping={};inside=False;header=None
    with gzip.open(path,'rt',encoding='utf-8') as f:
        for line in f:
            if line.startswith('ID\t'):header=line.rstrip().split('\t');continue
            if header is None or line.startswith(('!','#','^')):continue
            v=line.rstrip('\n').split('\t')
            if len(v)!=len(header):continue
            gene=mapping_row(v,header)
            if gene:mapping[v[0]]=gene
    assert 'SIGLEC1' in mapping.values()
    return mapping

def matrix_read(path,mapping):
    md={};data=[];inside=False;columns=None
    with gzip.open(path,'rt',encoding='utf-8') as f:
        for line in f:
            v=next(csv.reader([line],delimiter='\t'))
            if not v:continue
            if v[0].startswith('!Sample_'):md.setdefault(v[0][8:],[]).append(v[1:])
            if v[0]=='!series_matrix_table_begin':inside=True;continue
            if v[0]=='!series_matrix_table_end':break
            if inside:
                if columns is None:columns=v[1:];continue
                if v[0] in mapping:data.append([v[0],mapping[v[0]]]+list(map(float,v[1:])))
    d=pd.DataFrame(data,columns=['probe','gene']+columns)
    z=d.groupby('gene')[columns].median().T.reindex(columns=GENES)
    return z,md,d

def clariom():
    mapping=fetch_clariom_mapping()
    z,md,probes=matrix_read(IN/'GSE212067_series_matrix.txt.gz',mapping)
    ev.save(probes,'GSE212067_all_unique_target_probes')
    m=pd.DataFrame({'sample':md['geo_accession'][0],'title':md['title'][0]})
    m['patient']=m.title.str.extract(r'Patient (\d+)')[0]
    m['condition']=np.where(m.title.str.contains('Initial'),'primary','recurrent')
    for values in md.get('characteristics_ch1',[]):
        if values[0].startswith('age at diagnosis:'):m['age']=pd.to_numeric(pd.Series(values).str.split(': ',n=1).str[1],errors='coerce')
    assert m.age.ge(18).all();d=m.merge(z,left_on='sample',right_index=True,validate='one_to_one');d=ev.outcomes(d)
    ev.save(d,'GSE212067_specimen_values')
    ev.paired(d,'GSE212067','tissue','deposited normalized signal')

def u133plus():
    mapping=annot_mapping(IN/'GPL570.annot.gz')
    z,md,probes=matrix_read(IN/'GSE7696_series_matrix.txt.gz',mapping)
    ev.save(probes,'GSE7696_all_unique_target_probes')
    m=pd.DataFrame({'sample':md['geo_accession'][0],'title':md['title'][0]})
    for values in md.get('characteristics_ch1',[]):
        if not values[0]:continue
        key=values[0].split(':',1)[0];m[key]=[v.split(':',1)[1].strip() if ':' in v else '' for v in values]
    m['condition']=m['disease status'].map({'GBM':'primary','recurrent GBM':'recurrent'})
    d=m[m.condition.notna()].merge(z,left_on='sample',right_index=True,validate='one_to_one')
    if 'age' in d:d=d[pd.to_numeric(d.age,errors='coerce').ge(18)]
    # Average repeat specimens within the same author-provided patient/stage;
    # never infer pairing from similar sample titles.
    ev.save(d,'GSE7696_source_manifest')
    d=d.groupby(['patient','condition'],as_index=False)[GENES].mean()
    # Source study's processed RMA values are log-like; do not log-transform twice.
    d=ev.outcomes(d);ev.save(d,'GSE7696_specimen_values')
    ev.paired(d,'GSE7696','tissue','deposited normalized expression')

def phillips():
    mapping={};sample_platform={};data={};pid=None;sid=None;table=None;header=None
    with gzip.open(IN/'GSE4271_family.soft.gz','rt',encoding='utf-8') as f:
        for line in f:
            line=line.rstrip('\n')
            if line.startswith('^PLATFORM = '):pid=line.split(' = ')[1];mapping[pid]={}
            elif line.startswith('^SAMPLE = '):sid=line.split(' = ')[1];data[sid]={}
            elif line.startswith('!Sample_platform_id = '):sample_platform[sid]=line.split(' = ')[1]
            elif line.startswith('!platform_table_begin'):table='platform';header=None
            elif line.startswith('!sample_table_begin'):table='sample';header=None
            elif line.startswith(('!platform_table_end','!sample_table_end')):table=None
            elif table:
                fields=line.split('\t')
                if header is None:header=fields;continue
                if table=='platform':
                    g=mapping_row(fields,header)
                    if g:mapping[pid][fields[0]]=g
                elif fields[0] in mapping[sample_platform[sid]]:data[sid][fields[0]]=float(fields[1])
    md=json.loads((IN/'GSE4271_metadata.json').read_text())['samples'];rows=[]
    for s in md:
        platform=s['platform_id'][0]
        # Prefer the original U133A platform irrespective of expression outcomes.
        if platform!='GPL96':continue
        chars=dict(v.split(': ',1) for v in s['characteristics_ch1'] if ': ' in v)
        vals={g:[] for g in GENES}
        for probe,value in data[s['accession']].items():
            assert value>=0
            vals[mapping[platform][probe]].append(np.log2(value+1))
        if not vals['SIGLEC1']:raise ValueError('GPL96 lacks SIGLEC1')
        rows.append(dict(sample=s['accession'],title=s['title'][0],platform=platform,condition=chars['specimen type'],matching=chars.get('matching sample',''),grade=chars['who grade'],age=pd.to_numeric(chars['age'],errors='coerce'),**{g:np.median(v) if v else np.nan for g,v in vals.items()}))
    d=pd.DataFrame(rows);paired=[]
    for _,p in d[d.condition.eq('primary') & d.grade.str.startswith('IV') & d.age.ge(18)].iterrows():
        r=d[d.title.eq(p.matching)&d.condition.eq('recurrent')&d.grade.str.startswith('IV')]
        if len(r)!=1:continue
        r=r.iloc[0];assert r.matching==p.title
        for v in [p,r]:q=v.to_dict();q['patient']=p.title;paired.append(q)
    a=ev.outcomes(pd.DataFrame(paired));ev.save(d,'GSE4271_source_manifest');ev.save(a,'GSE4271_specimen_values')
    ev.paired(a,'GSE4271','tissue','log2(MAS5 signal + 1)')
    (ROOT/'GSE4271_probe_mapping.json').write_text(json.dumps(mapping,indent=2),encoding='utf-8')

if __name__=='__main__':
    status=[]
    for acc,fn in [('GSE212067',clariom),('GSE7696',u133plus),('GSE4271',phillips)]:
        try:fn();status.append(dict(cohort=acc,status='analyzed'))
        except Exception as e:status.append(dict(cohort=acc,status='error',error=str(e)));print(acc,type(e).__name__,str(e),flush=True)
    if ev.tests:
        d=pd.DataFrame(ev.tests);d['q_within_five']=d.groupby('cohort').p.transform(lambda v:ev.multipletests(v.fillna(1),method='fdr_bh')[1]);d.loc[d.p.isna(),'q_within_five']=np.nan;ev.save(d,'array_fixed_outcome_tests')
        print(d[d.metric.eq('SIGLEC1')].to_string(index=False),flush=True)
    (ROOT/'array_analysis_audit.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
