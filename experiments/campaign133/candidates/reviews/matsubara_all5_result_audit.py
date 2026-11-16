"""Bounded existing-state audit. stdlib+h5dump only; no Julia/new Phi/SCBA."""
from pathlib import Path
import resource,hashlib,tomllib,json,math,struct,tempfile,subprocess,re
resource.setrlimit(resource.RLIMIT_AS,(128*1024**2,128*1024**2))
p=Path('.build/campaign133/adjacent/matsubara-all5');o=p/'attempt01';read=lambda n:tomllib.loads((o/n).read_text());s=read('summary.toml');cfg=read('input.toml');ident=read('identity.toml');binding=read('source-binding.toml');ctrl=read('controls.toml');r={'summary_sha256':hashlib.sha256((o/'summary.toml').read_bytes()).hexdigest(),'checks':{}};c=r['checks'];near=lambda a,b:math.isclose(a,b,rel_tol=4e-12,abs_tol=1e-18)
frozen=tomllib.loads((p/'runtime-input.toml').read_text());c['runtime_only_auth_delta']={k for k in cfg.keys()|frozen.keys() if cfg.get(k)!=frozen.get(k)}=={'runtime_authorized'};c['input_sha']=hashlib.sha256((o/'input.toml').read_bytes()).hexdigest()==ident['input_sha256'];c['copied_source_hashes']=hashlib.sha256((o/'matsubara_all5.jl').read_bytes()).hexdigest()==cfg['text_sha256']['prototype'] and hashlib.sha256((o/'run_matsubara.jl').read_bytes()).hexdigest()==cfg['adapter_sha256'];c['source_frame']=ident['definition_sha256']==cfg['definition_sha256'] and ident['source_file_sha256']==cfg['source_file_sha256'] and binding['mu_native']==cfg['mu_native'] and binding['E0_eV']==.1 and binding['spin']==2
c['execution_counts']=s['execution']=='complete_bounded_diagnostic' and s['inner_requested']==s['inner_completed']==32 and s['completed_N']==[128,256] and s['status_N128']==s['status_N256']=='iteration_budget';c['scope_flags']=not any(s[k] for k in ('scientific_acceptance','native_pass','source_Gamma_used','density_normalized'))
rows={N:[read(f'N{N}-i{i:02}.toml') for i in range(1,17)] for N in (128,256)}
c['32_map_records']=all(len(rows[N])==16 and all(x['iteration']==i for i,x in enumerate(rows[N],1)) for N in rows);c['last_summary_matches']=s['last_record']==rows[256][-1];c['zero_iterations_11_to_16']=all(x['raw_sigma_absolute']==x['raw_sigma_relative']==0 for N in rows for x in rows[N][10:]);c['all_map_causality_and_Dyson']=all(x['sigma_causality_max']<0 and x['G_causality_max']<0 and x['dyson_backward']<1e-14 for N in rows for x in rows[N]);r['transient_pauli_negative_records']=[{'N':N,'iteration':x['iteration'],'pauli_min':x['pauli_min']} for N in rows for x in rows[N] if x['pauli_min']<0];c['native_all5_tiny_contractions']=set(ctrl['opcheck'])=={'LO','acoustic','impurity','IFR','alloy'} and max(ctrl['opcheck'].values())<2e-12
c['conditional_CP_ratio']=near(ctrl['cp']['L']/(math.pi/binding['beta_native'])**2,ctrl['cp']['conditional_L_over_omega0_squared'])
# Recompute independent scalar-oracle norm statements from preserved raw/exact complex values.
for row in ctrl['pole']:
 raw=complex(row['raw']['real'],row['raw']['imag']);exact=complex(row['exact']['real'],row['exact']['imag']);assert near(abs(raw-exact),row['absolute']) and near(abs(raw-exact)/abs(exact),row['relative'])
c['pole_refinement']=all(ctrl['pole'][i+1]['absolute']<ctrl['pole'][i]['absolute'] for i in (0,2,4));c['wired_refinement']=ctrl['wired'][1]['LO_absolute']<ctrl['wired'][0]['LO_absolute'] and ctrl['wired'][1]['density_absolute']<ctrl['wired'][0]['density_absolute']
requests=[read(f'request-N{N}.toml') for N in (128,256)];c['request_bounds']=requests[0]['total_requested']==16 and requests[0]['completed']==0 and requests[1]['total_requested']==32 and requests[1]['completed']==16
norm=lambda xs:math.sqrt(math.fsum(abs(z)**2 for z in xs));final={};r['independent_raw']={}
with tempfile.TemporaryDirectory(prefix='qcl-mats-result-') as td:
 def binary(f,name):
  dest=Path(td)/'dataset.bin';subprocess.run(['h5dump','-d',name,'-b','LE','-o',str(dest),str(f)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE);return dest.read_bytes()
 def doubles(f,name):
  buf=binary(f,name);return struct.unpack('<'+'d'*(len(buf)//8),buf)
 def field(f,name):
  a=doubles(f,name+'_real');b=doubles(f,name+'_imag');return [complex(x,y) for x,y in zip(a,b)]
 wk=doubles(Path('.build/campaign133/equilibrium-gamma-stableA-cont64/attempt01/spectral-selected-064.h5'),'weight_k');r['weight_k_binary_sha256']=hashlib.sha256(struct.pack('<11d',*wk)).hexdigest()
 for N in (128,256):
  full=o/f'N{N}-latest.h5';Sigma=field(full,'Sigma');Phi=field(full,'Phi');G=field(full,'G');rho=field(full,'rho');omega=doubles(full,'omega_native')
  exact=all(binary(full,'Sigma_'+part)==binary(full,'Phi_'+part) for part in ('real','imag'))
  number=2*math.fsum(wk[k]*rho[(a*5+a)*11+k].real for k in range(11) for a in range(5));hermdef=max(abs(rho[(b*5+a)*11+k]-rho[(a*5+b)*11+k].conjugate()) for k in range(11) for a in range(5) for b in range(5))
  # Independent floating complex LDL witnesses for rho and I-rho (not exact certificates).
  minpivot=1.;maxherm=0.
  for k in range(11):
   R=[[rho[(b*5+a)*11+k] for b in range(5)] for a in range(5)]
   for complement in (False,True):
    A=[[(1. if a==b else 0.)-R[a][b] if complement else R[a][b] for b in range(5)] for a in range(5)];L=[[0j]*5 for _ in range(5)];d=[]
    for j in range(5):
     pivot=A[j][j]-sum(L[j][h]*d[h]*L[j][h].conjugate() for h in range(j));assert abs(pivot.imag)<1e-14 and pivot.real>0
     d.append(pivot.real);minpivot=min(minpivot,pivot.real);L[j][j]=1
     for i in range(j+1,5):L[i][j]=(A[i][j]-sum(L[i][h]*d[h]*L[j][h].conjugate() for h in range(j)))/d[j]
  trajectory=[];previous=field(o/f'N{N}-i01.h5','Sigma');assert all(z==0 for z in previous)
  for i in range(1,16):
   nxt=field(o/f'N{N}-i{i+1:02}.h5','Sigma');absolute=norm([y-x for x,y in zip(previous,nxt)]);scale=norm(nxt);relative=absolute/scale if scale else 0.;rec=rows[N][i-1];assert near(absolute,rec['raw_sigma_absolute']) and near(relative,rec['raw_sigma_relative']),(N,i,absolute,rec)
   trajectory.append(relative);previous=nxt
  c[f'N{N}_latest_binary_Sigma_Phi_equal']=exact;c[f'N{N}_number_recomputed']=near(number,rows[N][-1]['number']);c[f'N{N}_Hermitian_rho']=hermdef<1e-14;c[f'N{N}_omega_grid']=all(near(w,(2*n+1)*math.pi/binding['beta_native']) for n,w in enumerate(omega))
  r['independent_raw'][str(N)]={'Sigma_Phi_bitwise_equal':exact,'raw_difference_norm':norm([x-y for x,y in zip(Sigma,Phi)]),'number':number,'rho_hermiticity_max':hermdef,'rho_and_IminusRho_floating_LDL_minpivot':minpivot,'successorSigma_15_raw_residual_replays':trajectory}
  final[N]={'G':G,'Sigma':Sigma,'rho':rho,'number':number}
 # Independent transient rho negative Rayleigh witness, tiny5x5 algebra only.
 rho2=field(o/'N128-i02.h5','rho');witness=None
 for kk in range(11):
  RR=[[rho2[(b*5+a)*11+kk] for b in range(5)] for a in range(5)];LL=[[0j]*5 for _ in range(5)];dd=[]
  for j in range(5):
   pivot=RR[j][j]-sum(LL[j][h]*dd[h]*LL[j][h].conjugate() for h in range(j));dd.append(pivot.real);LL[j][j]=1
   for i in range(j+1,5):LL[i][j]=(RR[i][j]-sum(LL[i][h]*dd[h]*LL[j][h].conjugate() for h in range(j)))/dd[j]
  if min(dd)<-1e-12:
   j=min(range(5),key=lambda z:dd[z]);v=[0j]*5
   for i in range(4,-1,-1):v[i]=(1. if i==j else 0.)-sum(LL[k][i].conjugate()*v[k] for k in range(i+1,5))
   den=sum(abs(x)**2 for x in v);q=sum(v[a].conjugate()*RR[a][b]*v[b] for a in range(5) for b in range(5))/den
   witness={'k_index':kk+1,'negative_LDL_pivot':dd[j],'normalized_Rayleigh_real':q.real,'Rayleigh_imag':q.imag,'vector':[[x.real,x.imag] for x in v]};assert q.real<-1e-12;break
 r['transient_independent_rho_witness']=witness;c['transient_negative_rho_raw_confirmed']=witness is not None
 # HDF flattened byte order is b,a,k,n. Compare actual first32 frequencies.
 def low(f,N):return [f[((b*5+a)*11+k)*N+n] for b in range(5) for a in range(5) for k in range(11) for n in range(32)]
 comparison={}
 for name,stem in [('G','G'),('Sigma','S')]:
  a=low(final[128][name],128);b=low(final[256][name],256);absolute=norm([x-y for x,y in zip(a,b)]);relative=absolute/norm(b);comparison[stem+'_absolute']=absolute;comparison[stem+'_relative']=relative
 comparison['rho_absolute']=norm([x-y for x,y in zip(final[128]['rho'],final[256]['rho'])]);comparison['number_difference']=final[128]['number']-final[256]['number']
 r['independent_cutoff']=comparison;c['cutoff_comparison_from_raw']=all(near(v,s['cutoff_comparison'][k]) for k,v in comparison.items())
r['density']={'N128':final[128]['number'],'N256':final[256]['number'],'cutoff_diff_over_target':abs(comparison['number_difference'])/.045,'N128_target_relative':(final[128]['number']-.045)/.045,'N256_target_relative':(final[256]['number']-.045)/.045};log=(p/'actual-root.log').read_text();rss=int(re.search(r'Maximum resident set size \(kbytes\): (\d+)',log)[1]);r['cost']={'process_seconds':25.12,'script_seconds':s['elapsed_seconds_including_imports'],'MaxRSS_KiB':rss,'MaxRSS_bytes':rss*1024,'swap0':'Swaps: 0' in log,'exit0':'Exit status: 0' in log,'output_bytes':sum(f.stat().st_size for f in o.rglob('*') if f.is_file())};c['cost_scope']=r['cost']['swap0'] and r['cost']['exit0'] and r['cost']['output_bytes']<=cfg['data_max_bytes'] and r['cost']['MaxRSS_bytes']==s['maxrss_bytes'];r['all_pass']=all(c.values());assert r['all_pass'],{a:b for a,b in c.items() if not b};Path('.build/campaign133/reviews/matsubara-all5-result-audit.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
