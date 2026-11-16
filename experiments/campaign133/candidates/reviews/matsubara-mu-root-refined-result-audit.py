"""7-call/112-map metadata and binary Sigma/Phi terminal audit; no new solver."""
import resource
resource.setrlimit(resource.RLIMIT_AS,(128*1024**2,128*1024**2))
from pathlib import Path
import hashlib,json,tomllib,math,subprocess,tempfile,array,sys,re
campaign=Path(__file__).resolve().parents[1];base=campaign/'adjacent/matsubara-mu-root-refined';p=base/'attempt02'
load=lambda q:tomllib.loads(q.read_text())
def sha(q):
 h=hashlib.sha256()
 with q.open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
s=load(p/'summary.toml');c=load(p/'input.toml')
assert sha(p/'summary.toml')=='e041cca82a7e9c08aa1a823fe24c09f906f87249d634c87cd4d11ab86f9e8bf3'
assert s['outcome']=='qualified_fixedcharge_M3_variant'
assert s['maps_completed_callbacks']==s['reserved_maps']==112<=s['structural_map_bound']==144
assert s['outer_calls']==6 and s['final_calls']==1 and len(s['calls'])==7
assert s['native_pass'] is s['scientific_acceptance'] is s['infinite_N_certified'] is False
assert s['BLAS_threads_final']==1
assert sha(p/'input.toml')=='0057daf02b9c5fd80f0384350d6c18de35ac07795ed0369dc9045175459a8191'
assert sha(p/'run_mu_root.jl')==c['wrapper_sha256']=='cb1d60a42c0be18761d5fa8e69846c1e5c8a81a08eb8b5c95ebc22483c31f144'
copies={'source_adapter':'source_adapter-run_matsubara.jl','source_module':'source_module-matsubara_all5.jl',
 'refined_module':'refined_module-matsubara_all5_refined.jl','M3_module':'M3_module-moment_tail.jl',
 'M3_review':'M3_review-matsubara-moments-result.md','M3_qualification':'M3_qualification-summary.toml',
 'M3_qualification_input':'M3_qualification_input-input.toml','source_input':'source_input-runtime-input.toml'}
for k,n in copies.items():assert sha(p/n)==c[k+'_sha256']
calls=[];last3=[]
for j,call in enumerate(s['calls'],1):
 N=512 if j<=6 else 1024;d=p/f'call-{j}-N{N}';hist=call['history']
 assert call['N']==N and call['S0']=='zero' and call['alpha']==1 and call['status']=='iteration_budget' and call['inner_qualified']
 assert len(hist)==16 and call['reserved_maps']==16
 assert len(list(d.glob('i*.toml')))==len(list(d.glob('i*.h5')))==16
 density=load(d/'density.toml');req=load(d/'request.toml');qual=load(d/'density-qualification.toml')
 assert req['mu_native']==density['mu_native']==call['mu_native'] and req['status']=='started_reserved'
 assert qual['candidate_accepted'] is True and density['candidate_accepted'] is False
 assert density['N']==N and density['eta']<1 and density['linear_residual_relative']<=1e-12
 for i,row in enumerate(hist,1):
  saved=load(d/f'i{i:02d}.toml');assert saved['raw']==row and saved['mu_native']==call['mu_native'] and saved['N']==N
  assert row['iteration']==i and all(math.isfinite(v) for v in row.values())
  assert saved['weighted_raw']['relative']==saved['weighted_raw']['absolute']/saved['weighted_raw']['scale']
  assert math.isclose(row['raw_sigma_relative'],row['raw_sigma_absolute']/saved['raw_unweighted_scale'],rel_tol=2e-14)
 for row in hist[-3:]:
  assert row['raw_sigma_relative']<=1e-9 and row['raw_sigma_absolute']<=1e-12;last3.append(row)
 calls.append({'call':j,'N':N,'mu':call['mu_native'],'number_M3':density['number'],
   'number_raw':density['raw_number'],'target_relative':density['target_relative'],
   'density_path':str(d/'density.toml')})
 assert math.isclose(density['target_relative'],(density['number']-c['target_native'])/c['target_native'],abs_tol=1e-20)
# Independently replay bounded sign-bracket/secant sequence of the six outer calls.
a,b=calls[:2];assert a['mu']==c['source_mu_native']-.001 and b['mu']==c['source_mu_native']+.001
assert a['target_relative']*b['target_relative']<0
for row in calls[2:6]:
 mu=b['mu']-b['target_relative']*(b['mu']-a['mu'])/(b['target_relative']-a['target_relative'])
 if not math.isfinite(mu) or not a['mu']<mu<b['mu']:mu=(a['mu']+b['mu'])/2
 assert mu==row['mu']
 if abs(row['target_relative'])<=c['charge_relative_tolerance']:break
 if math.copysign(1,row['target_relative'])==math.copysign(1,a['target_relative']):a=row
 else:b=row
assert calls[5]['mu']==calls[6]['mu']==s['selected']['mu_native']==s['witness']['mu_native']
assert load(p/'call-6-N512/density.toml')==s['selected'] and load(p/'call-7-N1024/density.toml')==s['witness']
assert abs(s['selected']['target_relative'])<=1e-9 and abs(s['witness']['target_relative'])<=1e-8
cutoff=(s['witness']['number']-s['selected']['number'])/c['target_native'];assert abs(cutoff)<=1e-8
binding=load(p/'source-binding.toml');beta=binding['beta_native'];binary=[]
def dump64(file,name,shape,tmp):
 raw=Path(tmp)/name
 header=subprocess.check_output(['h5dump','-H','-d','/'+name,str(file)],text=True,timeout=10)
 assert 'H5T_IEEE_F64LE' in header and tuple(map(int,re.search(r'SIMPLE \{ \( ([^)]*) \)',header).group(1).split(',')))==shape
 subprocess.run(['h5dump','-d','/'+name,'-o',str(raw),'-b','LE',str(file)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,timeout=15)
 count=math.prod(shape);assert raw.stat().st_size==8*count
 x=array.array('d')
 with raw.open('rb') as f:x.fromfile(f,count)
 if sys.byteorder!='little':x.byteswap()
 assert all(math.isfinite(v) for v in x)
 return x,sha(raw)
for j,N in [(6,512),(7,1024)]:
 d=p/f'call-{j}-N{N}';file=d/'latest.h5'
 with tempfile.TemporaryDirectory(prefix='mats-raw-audit-') as tmp:
  ar={};hashes={}
  for name in ['Sigma_real','Sigma_imag','Phi_real','Phi_imag']:
   ar[name],hashes[name]=dump64(file,name,(5,5,11,N),tmp)
  w,_=dump64(file,'omega_native',(N,),tmp)
  assert all(math.isclose(w[n],(2*n+1)*math.pi/beta,rel_tol=3e-16) for n in range(N))
  raw2=math.fsum((a-b)**2 for r,t in [('Sigma_real','Phi_real'),('Sigma_imag','Phi_imag')] for a,b in zip(ar[r],ar[t]))
  scale2=math.fsum(v*v for key in ('Phi_real','Phi_imag') for v in ar[key]);assert scale2>0
  absnorm=math.sqrt(raw2);relative=absnorm/math.sqrt(scale2);stored=s['calls'][j-1]['history'][-1]
  assert absnorm==stored['raw_sigma_absolute'] and relative==stored['raw_sigma_relative']
  binary.append({'call':j,'N':N,'latest_h5_sha256':sha(file),'raw_absolute':absnorm,'raw_relative':relative,
    'Phi_scale_independent':math.sqrt(scale2),'real_bitwise_equal':hashes['Sigma_real']==hashes['Phi_real'],
    'imag_bitwise_equal':hashes['Sigma_imag']==hashes['Phi_imag'],'datasets_binary_sha256':hashes,
    'logical_shape':[N,11,5,5],'HDF5_physical_shape':[5,5,11,N],'all_finite':True})
sizes=[q.stat().st_size for q in p.rglob('*') if q.is_file()]
assert sum(sizes)<=c['data_max_bytes'] and max(sizes)<=c['single_file_max_bytes']
result={'status':'all_metadata_and_two_binary_raw_roots_pass','summary_sha256':sha(p/'summary.toml'),
 'input_sha256':sha(p/'input.toml'),'wrapper_sha256':sha(p/'run_mu_root.jl'),'maps_completed':112,
 'outer_calls':6,'final_calls':1,'last3_max_relative':max(x['raw_sigma_relative'] for x in last3),
 'last3_max_absolute':max(x['raw_sigma_absolute'] for x in last3),'calls':calls,'binary_terminal_roots':binary,
 'mu_native':s['selected']['mu_native'],'mu_eV':s['selected']['mu_native']*binding['E0_eV'],
 'target_native':c['target_native'],'outer_target_relative':s['selected']['target_relative'],
 'witness_target_relative':s['witness']['target_relative'],'crossN_over_target':cutoff,
 'raw_N512_over_target':s['selected']['raw_number']/c['target_native']-1,
 'raw_N1024_over_target':s['witness']['raw_number']/c['target_native']-1,
 'data_bytes':sum(sizes),'no_Julia_newsolver_or_Gamma_source_read':True,
 'M3_value_Pauli_independent_audit':'separate b0_review artifact; not duplicated here'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('calls',)},indent=2))
