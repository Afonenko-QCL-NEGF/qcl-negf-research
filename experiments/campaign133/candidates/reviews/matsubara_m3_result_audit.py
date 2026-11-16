"""Bounded stdlib audit of saved moment matrices and literal K arithmetic. No Julia/physical Phi/Dyson."""
from pathlib import Path
import hashlib,json,math,re,resource,signal,struct,subprocess,time,tomllib
resource.setrlimit(resource.RLIMIT_AS,(128*1024**2,128*1024**2))
signal.alarm(30)
START=time.monotonic();BASE=Path(__file__).resolve().parents[1]
P=BASE/'adjacent/matsubara-moments';A=P/'attempt02';OUT=Path(__file__).resolve().parent
load=lambda p:tomllib.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(65536):h.update(b)
 return h.hexdigest()
def array_sha(v,complex_values=False):
 h=hashlib.sha256()
 for x in v:h.update(struct.pack('<dd',x.real,x.imag) if complex_values else struct.pack('<d',x))
 return h.hexdigest()
s=load(A/'summary.toml');c=load(A/'input.toml');frozen=load(P/'input.toml')
assert sha(A/'summary.toml')=='cb2dcc61903627ecd822f1152ba3d04c1f9679ed6c94344a28e6c95e0774ecc5'
assert c['runtime_authorized'] and not frozen['runtime_authorized']
assert dict(c,runtime_authorized=False)==frozen
assert s['execution']=='completed_postprocess_only' and s['inner_maps']==s['mu_fits']==0
assert s['BLAS_threads']==s['BLAS_threads_final']==s['Julia_threads']==1
for name,key in [('run_moments.jl','adapter'),('moment_tail.jl','moment_module')]:assert sha(A/name)==c['text_sha256'][key]
for key,h in c['text_sha256'].items():assert sha(P/c[key])==h
source=(P/c['actual_directory']).resolve()
for name,h in c['source_sha256'].items():assert sha(source/name)==h
old=load(source/'input.toml');packet=load(source/'source-binding.toml')
source_model=(source.parent/old['source_directory']).resolve()
assert sha(source_model/'model.h5')==old['source_file_sha256']['model.h5']
assert sha(source_model/'spectral-selected-064.h5')==old['source_file_sha256']['spectral-selected-064.h5']
commands=[]
def get(path,name):
 cmd=['h5dump','-y','-w0','-m','%.17g','-d',name,str(path)];commands.append(cmd)
 text=subprocess.run(cmd,check=True,capture_output=True,text=True,timeout=5).stdout
 data=re.findall(r'DATA\s*\{([^}]+)\}',text,re.S);assert len(data)==1
 result=[float(x) for x in data[0].replace(',',' ').split()];assert all(map(math.isfinite,result))
 return result
def cx(path,name,count=275):
 r,q=get(path,name+'_real'),get(path,name+'_imag');assert len(r)==len(q)==count
 return [complex(a,b) for a,b in zip(r,q)]
def norm(x):return math.sqrt(math.fsum(abs(z)**2 for z in x))
def sub(x,y):return [a-b for a,b in zip(x,y)]
def rel(a,b):return 0. if a==b==0 else a/b if b else math.inf
def ix(k,a,b):return k+11*(a+5*b)
model=source_model/'model.h5';wk=get(model,'grid_wᵏ');assert len(wk)==11
assert array_sha(wk)==packet['model_proof']['grid_wᵏ']
H=cx(source_model/'spectral-selected-064.h5','h_H');assert array_sha(H,True)==packet['hH_sha256']
mu=c['mu_native'];Om=packet['Omega_native'];beta=packet['beta_native'];NB=1/math.expm1(beta*Om)
I=[complex(a==b) for b in range(5) for a in range(5) for k in range(11)]
h=[x-mu*y for x,y in zip(H,I)]
def literal(K,q,X):
 # Julia K[m,mp,a,cc,d,b], X[mp,cc,d]. H5 C-flat order is Julia column-major.
 Y=[0j]*275
 for b in range(5):
  for a in range(5):
   for m in range(11):
    value=0j
    for mp in range(11):
     for cc in range(5):
      for d in range(5):
       index=m+11*(mp+11*(a+5*(cc+5*(d+5*b))))
       value+=wk[mp]*K[index]*X[ix(mp,cc,d)]
    Y[ix(m,a,b)]=q*value
 return Y
S0=[0j]*275;B0=[0j]*275;kernel_hashes={}
for family in ['acoustic','impurity','IFR','alloy','LO']:
 K=cx(model,f'kernels/{family}/K',75625);q=get(model,f'kernels/{family}/qK')[0]
 kh=array_sha(K,True);assert kh==packet['model_proof']['kernel_'+family]
 assert array_sha([q])==packet['model_proof']['qK_'+family]
 kernel_hashes[family]=kh
 ki=literal(K,q,I)
 if family=='LO':
  KLO,qLO=K,q
  khh=literal(K,q,[(2*NB+1)*x+Om*y for x,y in zip(h,I)])
  S0=[a+(2*NB+1)*b for a,b in zip(S0,ki)]
 else:
  khh=literal(K,q,h);S0=[a+b for a,b in zip(S0,ki)]
 B0=[a+b for a,b in zip(B0,khh)]
for k in range(11):
 for a in range(5):
  for b in range(5):
   B0[ix(k,a,b)]+=sum(h[ix(k,a,j)]*S0[ix(k,j,b)]+S0[ix(k,a,j)]*h[ix(k,j,b)] for j in range(5))
def sheet(X):return packet['spin']*math.fsum(wk[k]*math.fsum(X[ix(k,a,a)].real for a in range(5)) for k in range(11))
def wn(X):return math.sqrt(math.fsum(wk[k]*abs(X[ix(k,a,b)])**2 for b in range(5) for a in range(5) for k in range(11)))
def hermdef(X):return norm([X[ix(k,a,b)]-X[ix(k,b,a)].conjugate() for b in range(5) for a in range(5) for k in range(11)])
answer={'scope':'small saved densities and exact-source literal six-index K moment audit; no physical Phi/Dyson/new N',
 'source_hashes_verified':c['source_sha256'],'kernel_logical_hashes_verified':kernel_hashes,'N':{}}
fields={}
for N in [128,256]:
 rec=load(A/f'N{N}.toml');assert rec==s[f'N{N}']
 file=A/f'N{N}-densities.h5';assert get(file,'mu_native')==[mu]
 X={name:cx(file,name) for name in ['bare','classic_M1','oneshot','M3','feedback','interaction_S0','B0']}
 assert cx(source/f'N{N}-latest.h5','rho')==X['bare']
 dS0=rel(norm(sub(S0,X['interaction_S0'])),norm(S0));dB0=rel(norm(sub(B0,X['B0'])),norm(B0))
 assert dS0<2e-12 and dB0<2e-12
 t=rec['t']['t4'];a=2*Om*t;kr=literal(KLO,qLO,X['M3']);kb=literal(KLO,qLO,X['bare'])
 rhs=[x+t*y for x,y in zip(X['bare'],X['B0'])]
 residual=[x-y+a*z for x,y,z in zip(X['M3'],rhs,kr)]
 one=[x+t*(y-2*Om*z) for x,y,z in zip(X['bare'],X['B0'],kb)]
 item={'numbers':{name:sheet(X[name]) for name in ['bare','classic_M1','oneshot','M3']},
  'correction_Frobenius':norm(sub(X['M3'],X['bare'])),'feedback_Frobenius':norm(sub(X['M3'],X['oneshot'])),
  'stored_feedback_difference_Frobenius':norm(sub(X['feedback'],sub(X['M3'],X['oneshot']))),
  'oneshot_independent_literal_difference_Frobenius':norm(sub(one,X['oneshot'])),
  'S0_independent_literal_relative':dS0,'B0_independent_literal_relative':dB0,
  'linear_literal_residual_absolute':norm(residual),'linear_literal_residual_relative':norm(residual)/norm(rhs),
  'linear_literal_weighted_absolute':wn(residual),'linear_literal_weighted_relative':wn(residual)/wn(rhs),
  'M3_Hermiticity_defect':hermdef(X['M3']),'t4':t,'eta':rec['eta'],'condition':rec['linear_condition']}
 for name,key in [('bare','bare_number'),('classic_M1','classic_M1_number'),('oneshot','oneshot_number'),('M3','M3_number')]:assert math.isclose(item['numbers'][name],rec[key],abs_tol=3e-16)
 assert math.isclose(item['correction_Frobenius'],rec['correction_absolute'],rel_tol=1e-12)
 assert norm(residual)/norm(rhs)<1e-12
 assert item['stored_feedback_difference_Frobenius']==0
 answer['N'][str(N)]=item;fields[N]=X
answer['cutoff']={name:{'number_difference':sheet(fields[256][name])-sheet(fields[128][name]),
 'Frobenius_difference':norm(sub(fields[256][name],fields[128][name]))} for name in ['bare','classic_M1','oneshot','M3']}
answer['improvement_absolute_number_cutoff']=abs(answer['cutoff']['bare']['number_difference']/answer['cutoff']['M3']['number_difference'])
answer['M3_cutoff_relative_to_target']=abs(answer['cutoff']['M3']['number_difference'])/.045
answer['M3_N256_excess_over_target_relative']=(answer['N']['256']['numbers']['M3']-.045)/.045
v1=load(P/'attempt01/N128.toml');answer['V1_N128_reported_number_equal_V2']=all(v1[k]==s['N128'][k] for k in ['bare_number','classic_M1_number','oneshot_number','M3_number'])
answer['elapsed_seconds']=time.monotonic()-START;answer['peak_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
(OUT/'matsubara-m3-result-audit.json').write_text(json.dumps(answer,indent=2,allow_nan=False)+'\n')
(OUT/'matsubara-m3-result-h5commands.json').write_text(json.dumps(commands,indent=2)+'\n')
print(json.dumps({k:answer[k] for k in ['improvement_absolute_number_cutoff','M3_cutoff_relative_to_target','M3_N256_excess_over_target_relative','V1_N128_reported_number_equal_V2','elapsed_seconds','peak_KiB']},indent=2))
