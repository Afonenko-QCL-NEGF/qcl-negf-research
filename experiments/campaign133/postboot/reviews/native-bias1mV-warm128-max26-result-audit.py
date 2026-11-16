from pathlib import Path
import resource,hashlib,json,tomllib,csv,math,re,collections
resource.setrlimit(resource.RLIMIT_AS,(128*1024**2,128*1024**2))
r=Path(__file__).resolve().parents[1];p=r/'native-bias1mV-warm128-max26';a=p/'attempt01'
load=lambda x:tomllib.loads(x.read_text())
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
close=lambda x,y:math.isclose(x,y,rel_tol=2e-11,abs_tol=2e-13)
s=load(a/'summary.toml');c=load(a/'input.toml');e=json.loads((p/'actual-execution.json').read_text())
assert sha(a/'input.toml')=='c62fde3d99f3963d54e36bc6b8caeb9c056320983246e6add15383e749261919'
assert e['exit']==0 and s['execution']=='complete' and not s['scientific_acceptance'] and len(s['cases'])==2
assert c['runtime_authorized'] is True
out={'status':'all256_small_rows_replayed_both_unresolved','summary_sha256':sha(a/'summary.toml'),'source_sha256':sha(p/'run_native.jl'),'loader_sha256':sha(p/'warm_initial.jl'),'actual_input_sha256':sha(a/'input.toml'),'cases':[]}
for v in (0,1):
 d=a/f'case{v}';cs=load(d/'summary.toml');ci=load(d/'input.toml');id=load(d/'identity.toml');init=load(d/'warm-initialization.toml');bd=load(d/'boundary-replay.toml')
 assert cs==s['cases'][v] and cs['maps_completed']==128 and cs['status']=='max_iterations' and cs['quality']=='unresolved' and not cs['converged']
 assert ci['voltage_mV']==v and id['input_sha256']==sha(d/'input.toml')
 assert sha(d/'run_native.jl')==id['script_sha256']==out['source_sha256']
 assert sha(d/'warm_initial.jl')==id['warm_loader_sha256']==out['loader_sha256']
 assert sha(Path(id['project']).parent/'Manifest.toml')==id['manifest_sha256']
 assert init['parent_completed']==32 and not init['parent_converged'] and init['history_reset'] and init['consume_initial'] and not init['resume_iterations']
 parent=load(Path(init['parent'])/'summary.toml')['last_row']
 rows=[{k:float(x) for k,x in z.items()} for z in csv.DictReader((d/'history.csv').open())]
 assert len(rows)==128 and [x['ν'] for x in rows]==list(range(1,129))
 assert all(close(x,rows[-1][k]) for k,x in cs['last_row'].items())
 for row in rows:
  assert close(row['λ'],row['raw_charge']/row['target_charge']) and close(row['r_λ'],abs(row['λ']-1)) and close(row['normalized_charge'],row['target_charge'])
 assert bd['new_iteration']==1 and bd['parent_iteration']==32 and bd['green_replay_passed']
 for gf in bd['green_fields'].values():assert gf['bitwise_equal'] and all(gf[k]==0 for k in ('absolute','relative','weighted_absolute','weighted_relative'))
 for k,x in bd['fresh_equation_rows'].items():
  assert x['bitwise_equal'] and x['absolute_difference']==0 and x['fresh']==x['parent']==rows[0][k]==parent[k]
 expected=[1]+list(range(8,129,8));mi=sorted(int(x.stem.split('-')[1]) for x in d.glob('markers-*.toml'));bi=sorted(int(x.stem.split('-')[1]) for x in d.glob('balance-*.toml'))
 assert mi==bi==expected
 maxsig=0.
 for i in mi:
  row=rows[i-1];m=load(d/f'markers-{i:03d}.toml');b=load(d/f'balance-{i:03d}.toml')
  assert m['measured_iteration']==i and m['fresh_map_status']=='available'
  assert close(m['delta_energy_eV'],5e-5) and close(m['lo_shift_over_dE'],734) and close(m['field_shift_over_dE'],20*v)
  assert m['equilibrium_applicable']==(v==0)
  rr=max(x['residual_absolute']/(x['residual_scale']+1e-14) for x in m['channels'] if x['channel']!='total')
  assert close(rr,row['r_Σ']);maxsig=max(maxsig,abs(rr-row['r_Σ']))
  ii=[x for x in m['collisions'] if x['state_kind']=='mixed'];assert len(ii)==5
  signed=sum(x['particle_signed'] for x in ii)+b['plus']['particle_signed']+b['minus']['particle_signed']
  scale=sum(x['particle_absolute'] for x in ii)+b['plus']['particle_absolute']+b['minus']['particle_absolute']
  assert close(signed,b['particle_signed']) and close(scale,b['particle_scale']) and close(abs(signed)/scale,b['particle_relative'])
  assert close(row['J']*v*1e-3,b['Pfield_W_m2']) and close(abs(b['Psc_W_m2']+b['Pfield_W_m2']),b['power_balance_abs_W_m2'])
  assert close(b['power_balance_abs_W_m2']/b['power_scale_W_m2'],b['power_relative_diagnostic'])
 final=rows[-1];balance=load(d/'balance-128.toml');joint=[max(z['r_K'],z['r_Σ'],z['r_λ']) for z in rows]
 out['cases'].append({'V_mV_period':v,'maps':128,'novel_post_parent_states':127,'nominal_final_cumulative_index':159,'last_raw':final,'first_raw_matches_parent':True,'boundary_Green_bitwise_equal_reported':True,'marker_balance_iterations':mi,'rSigma_replay_maxabs':maxsig,'minimum_joint_raw':min(joint),'final_J_A_cm2_UNCONVERGED':final['J']/1e4,'particle_relative_diagnostic':balance['particle_relative'],'power_relative_diagnostic':balance['power_relative_diagnostic'],'elapsed_case_seconds':cs['elapsed_seconds'],'model_sha256_declared':cs['model_sha256'],'state_sha256_declared':cs['state_sha256']})
# Parse Julia's repr as metadata using regex only, never eval.
hdr=re.compile(r'^SolverEvent\(:(\w+), :(\w+), "([^"]*)", :(\w+), ([^,]+),')
metric=re.compile(r'SolverMetric\(:(\w+), ("(?:[^"\\]|\\.)*"|[^,]+), "([^"]*)"\)')
active={};case=-1;agg=collections.defaultdict(lambda:{'count':0,'wall':0.,'cpu_process':0.,'allocated_bytes_process':0.});progress=[]
for line in (p/'phase-events.log').open():
 h=hdr.match(line);assert h
 kind,stage,label,status,it=h.groups();it=int(it)
 if kind=='begin' and stage=='scba':case+=1;assert it==0 and ', 128,' in line
 if kind not in ('phase_begin','phase_end','progress'):continue
 vals={k:(json.loads(x) if x.startswith('"') else (x=='true' if x in ('true','false') else float(x))) for k,x,u in metric.findall(line)}
 if kind=='progress':
  assert 1<=it<=128;progress.append((case,it));assert vals['fixed_point_passed'] is False and vals['convergence_streak']==0
  continue
 assert vals['counter_scope']=='process' and 0<=it<=128
 key=(case,label,it)
 if kind=='phase_begin':assert not active;active[key]=vals['start_monotonic_ns'];continue
 assert status=='completed' and key in active and vals['start_monotonic_ns']==active.pop(key)
 assert close(vals['duration_seconds'],(vals['end_monotonic_ns']-vals['start_monotonic_ns'])/1e9)
 assert 0<=vals['duration_seconds']<=e['seconds'] and vals['cpu_seconds']>=0
 z=agg[label];z['count']+=1;z['wall']+=vals['duration_seconds'];z['cpu_process']+=vals['cpu_seconds'];z['allocated_bytes_process']+=vals['allocated_bytes']
assert case==1 and not active and progress==[(v,i) for v in (0,1) for i in range(1,129)]
wall=sum(x['wall'] for x in agg.values());cpu=sum(x['cpu_process'] for x in agg.values());stat=dict(z.split() for z in e['cgroup']['cpu.stat'].splitlines());cpuall=int(stat['usage_usec'])/1e6
assert wall<=e['seconds'] and cpu<=cpuall
out['phases']={'case_begins':2,'raw_progress_rows':len(progress),'balanced_begin_end':True,'counter_scope':'process','summed_phase_wall_seconds':wall,'summed_phase_cpu_seconds':cpu,'by_phase':dict(agg)}
out['resources']=e;out['resources']['average_cpu_cores']=cpuall/e['seconds'];out['resources']['fraction_of_26core_quota']=cpuall/e['seconds']/26
out['output_bytes']=s['bytes'];out['no_Julia_Pkg_newmaps_or_large_HDF_reads']=True
(r/'reviews/native-bias1mV-warm128-max26-result-audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'summary_sha256':out['summary_sha256'],'cases':[{k:x[k] for k in ('V_mV_period','maps','rSigma_replay_maxabs','minimum_joint_raw','particle_relative_diagnostic','power_relative_diagnostic')} for x in out['cases']],'phase_wall':wall,'phase_cpu':cpu,'CPUcores':cpuall/e['seconds']},indent=2))
