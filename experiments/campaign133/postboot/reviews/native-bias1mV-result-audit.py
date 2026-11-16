from pathlib import Path
import resource,hashlib,json,tomllib,csv,math
resource.setrlimit(resource.RLIMIT_AS,(128*1024**2,128*1024**2))
r=Path(__file__).resolve().parents[1];p=r/'native-bias1mV-postreboot';a=p/'attempt01'
load=lambda x:tomllib.loads(x.read_text())
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
close=lambda x,y:math.isclose(x,y,rel_tol=2e-11,abs_tol=2e-13)
s=load(a/'summary.toml');cfg=load(a/'input.toml');assert sha(a/'input.toml')=='34d184c9543a5cdd097efba89dc64dffbf265835dedb2029c55516ef891767d0'
assert s['execution']=='complete' and not s['scientific_acceptance'] and len(s['cases'])==2
assert cfg['runtime_authorized'] is True
out={'status':'small_metadata_replay_pass_no_large_state_reads','summary_sha256':sha(a/'summary.toml'),'input_sha256':sha(a/'input.toml'),'cases':[]}
for v in (0,1):
 d=a/f'case{v}';cs=load(d/'summary.toml');ci=load(d/'input.toml');ident=load(d/'identity.toml')
 assert cs==s['cases'][v] and cs['maps_completed']==32 and cs['status']=='max_iterations' and cs['quality']=='unresolved' and not cs['converged']
 assert ci['voltage_mV']==v and ci['runtime_authorized'] is True and ident['input_sha256']==sha(d/'input.toml')
 assert sha(d/'run_native.jl')==ident['script_sha256']==sha(p/'run_native.jl')
 assert sha(Path(ident['project']).parent/'Manifest.toml')==ident['manifest_sha256']
 assert sha(d/'model.h5')==cs['model_sha256'] # small7.5MiB only; no array decode
 rows=[{k:float(x) for k,x in row.items()} for row in csv.DictReader((d/'history.csv').open())]
 assert len(rows)==32 and [x['ν'] for x in rows]==list(range(1,33))
 assert all(close(x,rows[-1][k]) for k,x in cs['last_row'].items())
 maxsig=0.;maxbalance=0.
 for i,row in enumerate(rows,1):
  m=load(d/f'markers-{i:03d}.toml');b=load(d/f'balance-{i:03d}.toml')
  assert m['measured_iteration']==i and m['fresh_map_status']=='available'
  assert close(m['delta_energy_eV'],5e-5) and close(m['lo_shift_over_dE'],734) and close(m['field_shift_over_dE'],20*v)
  assert m['equilibrium_applicable']==(v==0) and m['equilibrium_status']==('available' if v==0 else 'finite_bias')
  assert close(row['λ'],row['raw_charge']/row['target_charge']) and close(row['r_λ'],abs(row['λ']-1))
  assert close(row['normalized_charge'],row['target_charge'])
  ch=[x for x in m['channels'] if x['channel']!='total']
  rr=max(x['residual_absolute']/(x['residual_scale']+1e-14) for x in ch)
  maxsig=max(maxsig,abs(rr-row['r_Σ']));assert close(rr,row['r_Σ'])
  internal=[x for x in m['collisions'] if x['state_kind']=='mixed'];assert len(internal)==5
  signed=sum(x['particle_signed'] for x in internal)+b['plus']['particle_signed']+b['minus']['particle_signed']
  scale=sum(x['particle_absolute'] for x in internal)+b['plus']['particle_absolute']+b['minus']['particle_absolute']
  assert close(signed,b['particle_signed']) and close(scale,b['particle_scale']) and close(abs(signed)/scale,b['particle_relative'])
  assert close(row['J']*v*1e-3,b['Pfield_W_m2']) and close(abs(b['Psc_W_m2']+b['Pfield_W_m2']),b['power_balance_abs_W_m2'])
  assert close(b['power_balance_abs_W_m2']/b['power_scale_W_m2'],b['power_relative_diagnostic'])
 last=rows[-1];m=load(d/'markers-032.toml');b=load(d/'balance-032.toml')
 out['cases'].append({'Vperiod_mV':v,'maps':32,'status':cs['status'],'quality':cs['quality'],'last_raw':last,'rSigma_marker_replay_maxabs':maxsig,'particle_relative_diagnostic':b['particle_relative'],'power_relative_diagnostic':b['power_relative_diagnostic'],'normalized_Gn_correction':m['relative_correction_Gn'],'sampled_linewidth_q10_q50_q90':[m['sampled_gamma_over_dE_q10'],m['sampled_gamma_over_dE_q50'],m['sampled_gamma_over_dE_q90']],'sampled_underresolved_fraction':m['sampled_spectral_weight_underresolved'],'final_J_A_cm2_UNCONVERGED':last['J']/1e4,'model_sha256_independent':sha(d/'model.h5'),'state_sha256_declared_not_rehashed':cs['state_sha256'],'elapsed_seconds':cs['elapsed_seconds']})
e=json.loads((p/'actual-execution.json').read_text());assert e['exit']==0
out['execution']=e;out['data_bytes']=s['bytes'];out['no_Julia_new_maps_or_large_state_arrays']=True
(r/'reviews/native-bias1mV-result-audit.json').write_text(json.dumps(out,indent=2,allow_nan=True)+'\n')
print(json.dumps({'status':out['status'],'cases':[{k:x[k] for k in ('Vperiod_mV','maps','quality','rSigma_marker_replay_maxabs','particle_relative_diagnostic','power_relative_diagnostic','final_J_A_cm2_UNCONVERGED')} for x in out['cases']],'summary_sha256':out['summary_sha256']} ,indent=2))
