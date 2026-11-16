"""Small certificate metadata plus two independent exact rational LDL proofs."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,tomllib,resource,math
resource.setrlimit(resource.RLIMIT_AS,(128*1024**2,128*1024**2))
campaign=Path(__file__).resolve().parents[1];p=campaign/'global-gamma-exact-pd'
load=lambda q:tomllib.loads(q.read_text());sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
s=load(p/'all-attempt01/summary.toml');c=load(p/'all-run.toml');root=json.loads((campaign/'stableA-cont64-source-digests.json').read_text())
assert s['mode']==c['mode']=='all' and s['execution']=='complete'
assert s['all_blocks_strict_PD_certified'] is True
assert s['results']['checked_blocks']==s['expected_all_blocks']==s['requested_blocks']==30704*11==337744
assert s['results']['counts']=={'pd':337744,'nonpositive':0,'nonhermitian':0,'nonfinite':0}
assert s['results']['checked_exact_divrem_components']==337744*60==20264640
assert not s['results']['failure_witnesses_first8'] and 'failure' not in s
assert s['script_sha256']==sha(p/'run_certificate.jl') and s['math_sha256']==sha(p/'exact_pd.jl')
assert s['input_sha256']==sha(p/'all-run.toml')
names={'gamma':'Gamma-last.h5','model':'model.h5','input':'input.toml','identity':'identity.toml','map':'map-064.toml','operator':'gamma_functions.jl','spectral':'stable_spectral.jl','ledger':'spectral-ledger-064.toml','summary':'summary.toml'}
for k,n in names.items():
    assert s['source_sha256'][k]==c['bindings'][k]['sha256']==root[n]
    if n.endswith(('.toml','.jl')):assert sha((p/c['bindings'][k]['path']).resolve())==root[n]
bench=load(p/'benchmark-attempt01/summary.toml')
assert bench['mode']=='benchmark5000' and bench['results']['checked_blocks']==5000 and not bench['all_blocks_strict_PD_certified']
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def mul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def conj(a):return (a[0],-a[1])
def scale(a,x):return (a[0]*x,a[1]*x)
zero=(F(0),F(0));one=(F(1),F(0)); proofs=[]
for w in s['results']['first_block']+[s['results']['smallest_physical_determinant_positive_block']]:
    assert w['status']=='pd' and w['stage']==5
    Q=[[(F(w['Gamma_real'][5*i+j]),F(w['Gamma_imag'][5*i+j])) for j in range(5)] for i in range(5)]
    assert all(Q[i][j]==conj(Q[j][i]) for i in range(5) for j in range(5))
    L=[[zero for _ in range(5)] for _ in range(5)];d=[];det=F(1)
    for j in range(5):
        pivot=Q[j][j]
        for k in range(j):pivot=sub(pivot,scale(mul(L[j][k],conj(L[j][k])),d[k]))
        assert pivot[1]==0 and pivot[0]>0;d.append(pivot[0]);L[j][j]=one
        det*=d[j];assert det==F(int(w['leading_minors_integer'][j]))*(F(2)**((j+1)*w['integer_scale_exponent']))
        for i in range(j+1,5):
            value=Q[i][j]
            for k in range(j):value=sub(value,scale(mul(L[i][k],conj(L[j][k])),d[k]))
            L[i][j]=scale(value,1/d[j])
    for i in range(5):
        for j in range(5):
            rebuilt=zero
            for k in range(5):rebuilt=add(rebuilt,scale(mul(L[i][k],conj(L[j][k])),d[k]))
            assert rebuilt==Q[i][j]
    proofs.append({'e':w['energy_index'],'k':w['k_index'],'exact_Hermitian':True,'positive_LDL_pivots':5,'all_five_published_minors_match':True,'exact_25entry_reconstruction':True,'determinant_log10_approx':math.log10(det.numerator)-math.log10(det.denominator)})
result={'status':'ALL_certificate_metadata_and_two_exact_witnesses_pass','summary_sha256':sha(p/'all-attempt01/summary.toml'),'input_sha256':s['input_sha256'],'math_sha256':s['math_sha256'],'script_sha256':s['script_sha256'],'source_Gamma_sha256':root['Gamma-last.h5'],'checked_blocks':337744,'component_divrem_checks':20264640,'proofs':proofs,'benchmark_excluded_from_ALL_count':True,'no_source_HDF5_or_Julia':True,'scope':'full-count fixed-code exact certificate; independent rational reconstruction only two saved matrices, not a second whole-source scan'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
