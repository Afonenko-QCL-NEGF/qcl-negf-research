"""Small stdlib TOML/script audit and selected H5dump only. No Julia/maps/full arrays."""
from pathlib import Path
import collections, csv, hashlib, json, re, struct, subprocess, tomllib

ROOT=Path(__file__).resolve().parents[1]
OUT=Path(__file__).resolve().parent
def load(p): return tomllib.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest() # scripts/TOML only
def bits(x): return struct.pack('>d',x).hex()
def bitwise_equal(a,b):
    if isinstance(a,float) and isinstance(b,float): return bits(a)==bits(b)
    if isinstance(a,dict) and isinstance(b,dict): return a.keys()==b.keys() and all(bitwise_equal(a[k],b[k]) for k in a)
    if isinstance(a,list) and isinstance(b,list): return len(a)==len(b) and all(bitwise_equal(x,y) for x,y in zip(a,b))
    return type(a)==type(b) and a==b
def data_dump(path,dataset,start=None):
    cmd=['h5dump','-y','-w0','-m','%.17g','-d',dataset]
    if start is not None: cmd+=['-s',','.join(map(str,start)),'-c','5,5,1,1']
    cmd+=[str(path)]
    text=subprocess.run(cmd,check=True,capture_output=True,text=True).stdout
    data=re.findall(r'DATA\s*\{([^}]+)\}',text,re.S)
    if len(data)!=1: raise RuntimeError('unexpected H5dump data frame')
    return cmd,text,[float(x) for x in data[0].replace(',',' ').split()]

result={'scope':'read-only small TOML/source + selected 5x5 Gamma H5dump; no full array equality asserted'}
runs=[]
for label,folder in [('control','equilibrium-gamma-stableA-cont64'),('projected','equilibrium-gamma-projected64')]:
    p=ROOT/folder/'attempt01';ident=load(p/'identity.toml');cfg=load(p/'input.toml');summary=load(p/'summary.toml')
    maps=[load(f) for f in sorted(p.glob('map-*.toml'))]
    attempts=[load(f) for f in sorted(p.glob('attempt-*.toml'))]
    projections=[load(f) for f in sorted(p.glob('projection-*.toml'))]
    sources={n:sha(p/n)==h for n,h in ident['source_sha256'].items()}
    assert all(sources.values()) and sha(p/'input.toml')==ident['input_sha256']
    assert maps[-1]==summary['reduced']['last_map']
    assert load(p/'spectral-ledger-064.toml')==maps[-1]['spectral_operator_ledger']
    counts=collections.Counter((a['reason'],a['kind']) for a in attempts)
    item={'input_sha256':sha(p/'input.toml'),'identity_sha256':sha(p/'identity.toml'),
          'summary_sha256':sha(p/'summary.toml'),'map64_sha256':sha(p/'map-064.toml'),
          'sources_match_identity':sources,'completed_map_records':len(maps),'map_calls':[m['call'] for m in maps],
          'attempt_count':len(attempts),'attempt_reason_kind_counts':{' | '.join(k):v for k,v in counts.items()},
          'final_residual':maps[-1]['residual'],'best':min((m['residual']['relative'],m['call']) for m in maps),
          'mu':maps[-1]['mu'],'summary_elapsed_seconds':summary['elapsed_seconds'],
          'precheck_actual_calls':summary['precheck_actual_calls'],'pi_calls':summary['pi_calls'],
          'native':summary['native'],'native_accounting':summary['native_accounting'],
          'trigger_count':sum(m['trigger'] for m in maps),
          'min_computed_input_eigenvalue':min(m['input_cone']['minimum_eigenvalue'] for m in maps),
          'min_computed_rawPhi_eigenvalue':min(m['map_cone']['minimum_eigenvalue'] for m in maps)}
    boundary=load(p/'boundary-replay.toml');assert boundary['passed'] and boundary['all_nine_scalars_bitwise_equal']
    assert boundary['raw_Phi_field_replay']['absolute']==0
    item['boundary']=boundary
    if projections:
        changed=[r for r in projections if r['changed_blocks']]
        item['projections']={'count':len(projections),'changed_count':len(changed),'changed_next_calls':[r['next_pi_call'] for r in changed],
             'changed_blocks':[r['changed_blocks'] for r in changed],
             'max_absolute':max(r['absolute'] for r in projections),'max_relative':max(r['relative'] for r in projections),
             'max_weighted_absolute':max(r['weighted_absolute'] for r in projections),
             'max_weighted_relative':max(r['weighted_relative'] for r in projections),
             'max_antiHerm':max(r['input_antihermitian_norm'] for r in projections),
             'max_eigen_backward_relative':max(r['maximum_eigen_backward_relative'] for r in projections),
             'max_eigen_orthogonality':max(r['maximum_eigenvector_orthogonality_defect'] for r in projections),
             'seconds':sum(r['seconds'] for r in projections),'first_changed':changed[0]}
        assert all(a['AA']['projection']==projections[a['AA']['projection']['projection_attempt']-1] for a in attempts if a.get('AA',{}).get('projection'))
        assert all(r['changed_blocks']+r['unchanged_exact_feasible_blocks']==30704*11 for r in projections)
        assert all(sum(r['output_positive_rank_counts'])==30704*11 for r in projections)
    with (OUT/f'projected64-{label}-maps.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['call','raw_relative','weighted_raw_relative','raw_absolute','weighted_raw_absolute','mu','trigger'])
        for m in maps:writer.writerow([m['call'],m['residual']['relative'],m['residual']['weighted_relative'],m['residual']['absolute'],m['residual']['weighted_absolute'],m['mu']['mu'],m['trigger']])
    result[label]=item;runs.append((p,cfg,maps,projections))

a,b=runs
result['input_deltas']={k:[a[1].get(k),b[1].get(k)] for k in sorted(set(a[1])|set(b[1])) if a[1].get(k)!=b[1].get(k)}
keys=['residual','AA_working_residual','mu','input_cone','map_cone','channel_Gamma_norms','spectral_operator_ledger']
result['first_four_scientific_scalar_fields_equal']={str(i+1):{k:bitwise_equal(a[2][i][k],b[2][i][k]) for k in keys} for i in range(4)}
assert all(all(checks.values()) for checks in result['first_four_scientific_scalar_fields_equal'].values())
for key in ['relative','weighted_relative']:
    result['improvement_'+key]=1-b[2][-1]['residual'][key]/a[2][-1]['residual'][key]

# HDF5 C order reverses the logical Julia (E,k,a,b) dimensions.
raw_dumps=[];checks=[];p=b[0]
for e,k in [(382,11),(5090,1),(30704,1)]:
    record={'energy_index':e,'k_index':k,'datasets':{}}
    for dataset in ['Gamma_real','Gamma_imag']:
        blocks=[]
        for name in ['Gamma-last.h5','Gamma-best.h5','reduced-terminal.h5']:
            cmd,text,values=data_dump(p/name,dataset,[0,0,k-1,e-1]);raw_dumps.extend([' '.join(cmd),text])
            assert len(values)==25;blocks.append([bits(x) for x in values])
        assert blocks[0]==blocks[1]==blocks[2]
        record['datasets'][dataset]={'selected_bitwise_equal_last_best_terminal':True,'last_bits':blocks[0]}
    checks.append(record)
cmd,text,values=data_dump(p/'reduced-terminal.h5','mu');raw_dumps.extend([' '.join(cmd),text])
assert len(values)==1 and bits(values[0])==bits(b[2][-1]['mu']['mu'])
result['selected_terminal_alignment']=checks;result['terminal_mu_bitwise_equal_map64']=True
(OUT/'projected64-h5-selected.txt').write_text('\n'.join(raw_dumps))
(OUT/'projected64-audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:result[k] for k in ['input_deltas','improvement_relative','improvement_weighted_relative','terminal_mu_bitwise_equal_map64']},indent=2))
