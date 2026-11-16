"""Independent small-artifact audit; no Julia, full Gamma/model, or full response arrays."""
from pathlib import Path
import resource, hashlib, tomllib, json, math, struct, subprocess, tempfile
resource.setrlimit(resource.RLIMIT_AS,(128*1024**2,128*1024**2))
p=Path('.build/campaign133/streaming-spectral');out=p/'qualified-attempt01'
s=tomllib.loads((out/'qualifier-summary.toml').read_text()); cfg=tomllib.loads((out/'input.toml').read_text())
r={'summary_sha256':hashlib.sha256((out/'qualifier-summary.toml').read_bytes()).hexdigest(),'checks':{},'selected_independent':{}}
c=r['checks'];c['execution']=s['execution']=='complete_operator_equivalence_only' and s['stage']=='complete';c['counters']=s['responses_started']==2 and s['old_response_completed'] and s['stream_response_completed'] and s['Phi_calls']==s['mu_fits']==s['mixing_steps']==0
c['immutable_flags']=all(s[k] for k in ('sources_unchanged_after','definitions_unchanged_after','Gamma_input_unchanged_after'));c['input_sha']=hashlib.sha256((out/'input.toml').read_bytes()).hexdigest()==s['input_sha256'];c['runtime_policy']=cfg['runtime_authorized'] and cfg['memory_max_bytes']==13019119616 and cfg['data_max_bytes']==32*1024**2
c['copied_code_hashes']=all(hashlib.sha256((out/n).read_bytes()).hexdigest()==h for n,h in s['code_sha256'].items())
root=json.loads(Path('.build/campaign133/stableA-cont64-source-digests.json').read_text());c['source_root_bindings']=all(root[n]==h==cfg['source_sha256'][n] for n,h in s['source_sha256'].items())
rows={f'{group}.{name}':v for group in ('all_entry_comparisons','selected_comparisons') for name,v in s[group].items()}
zeros=('absolute','weighted_absolute','relative','weighted_relative','maximum_block_absolute','maximum_block_relative','maximum_block_scaled')
c['all_15_comparisons_zero']=len(rows)==15 and all(all(v[k]==0 for k in zeros) and v['bitwise_equal'] for v in rows.values())
c['all_30_gates']=all(v['original_arithmetic_pass'] and v['planned_tight_pass'] for v in rows.values()) and s['original_arithmetic_pass'] and s['planned_tight_pass']
c['roundoff']=s['old_PI_roundoff']==s['stream_PI_roundoff'];c['memory_arithmetic']=s['memory_array_ledger']['paired_stream_sweep_6F']==6*135097600 and s['memory_array_ledger']['stream_own_sweep_4F']==4*135097600
with tempfile.TemporaryDirectory(prefix='qcl-stream-audit-') as td:
 def binary(file,name):
  dest=Path(td)/'dataset.bin'
  subprocess.run(['h5dump','-d',name,'-b','LE','-o',str(dest),str(file)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
  return dest.read_bytes()
 def doubles(file,name):
  buf=binary(file,name);return struct.unpack('<'+'d'*(len(buf)//8),buf)
 old=out/'spectral-selected-000.h5';new=out/'stream-selected.h5'
 chosen=s['selected_energy_indices'];ne=len(chosen);nk=11;nb=5
 wefull=doubles(old,'weight_energy');wk=doubles(old,'weight_k');we=[wefull[e-1] for e in chosen]
 for name in ('Gamma','Lambda_raw','Lambda_H_preassembly','GR','A_cov'):
  bufs=[];equal=True
  for suffix in ('real','imag'):
   x=binary(old,name+'_'+suffix);y=binary(new,name+'_'+suffix);equal &= x==y;bufs.append(struct.unpack('<'+'d'*(len(x)//8),x))
  re,im=bufs;ss=sw=0.0
  for b in range(nb):
   for a in range(nb):
    for k in range(nk):
     for e in range(ne):
      i=(((b*nb+a)*nk+k)*ne+e);v=re[i]**2+im[i]**2;ss+=v;sw+=we[e]*wk[k]*v
  v=s['selected_comparisons'][name];norm=math.sqrt(ss);weighted=math.sqrt(sw)
  ok=equal and math.isclose(norm,v['scale'],rel_tol=2e-13,abs_tol=0) and math.isclose(weighted,v['weighted_scale'],rel_tol=2e-13,abs_tol=0)
  r['selected_independent'][name]={'binary_equal':equal,'scale':norm,'weighted_scale':weighted,'norms_match_report':ok};assert ok,name
c['selected_5_fields_independent']=all(x['norms_match_report'] for x in r['selected_independent'].values())
r['all_pass']=all(c.values());r['full_array_scope']='5 full comparisons validated from fixed comparator code/report, not independent full-array reloading; only 5 selected complex fields independently reread'
r['scope']={'selected_E':s['selected_energy_indices'],'selected_blocks':len(s['selected_energy_indices'])*11,'responses':2,'Phi':0,'mu_fits':0,'maxrss_bytes_summary':s['maxrss_bytes'],'script_seconds':s['seconds'],'output_bytes':sum(f.stat().st_size for f in out.iterdir() if f.is_file())}
assert r['all_pass'],r
Path('.build/campaign133/reviews/streaming-qualifier-result-audit.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
