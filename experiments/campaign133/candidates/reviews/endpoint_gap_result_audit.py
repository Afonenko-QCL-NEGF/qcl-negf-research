"""Independent exact-rational reconstruction from small saved intervals/slices; no Julia/new G/Phi."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,tomllib,struct,math,tempfile,subprocess,resource,re
resource.setrlimit(resource.RLIMIT_AS,(128*1024**2,128*1024**2))
p=Path('.build/campaign133/adjacent/endpoint-gap-certificate');o=p/'attempt01';read=lambda n:tomllib.loads((o/n).read_text());s=read('summary.toml');cfg=read('input.toml');vc=read('V4-input.toml');vs=read('V4-summary.toml');r={'summary_sha256':hashlib.sha256((o/'summary.toml').read_bytes()).hexdigest(),'checks':{}};c=r['checks'];Q=F.from_float
frozen=tomllib.loads((p/'input.toml').read_text());c['runtime_only_auth_delta']={k for k in cfg.keys()|frozen.keys() if cfg.get(k)!=frozen.get(k)}=={'execution_authorized'} and cfg['execution_authorized'];c['hashes']=hashlib.sha256((o/'input.toml').read_bytes()).hexdigest()==s['input_sha256'] and hashlib.sha256((o/'endpoint_certificate.jl').read_bytes()).hexdigest()==s['script_sha256'] and hashlib.sha256((o/'V4-matsubara_oracle.jl').read_bytes()).hexdigest()==vc['oracle_script_sha256'];c['V4_actual_binding']=vs['status']=='completed' and s['V4_input_sha256']==hashlib.sha256((o/'V4-input.toml').read_bytes()).hexdigest()==cfg['V4_input_sha256'] and s['V4_summary_sha256']==hashlib.sha256((o/'V4-summary.toml').read_bytes()).hexdigest()==cfg['V4_summary_sha256'];c['complete_scope']=s['status']=='completed' and s['all_gaps_positive'] and not s['scientific_acceptance'] and s['Julia_threads']==s['BLAS_threads_after_setter']==1;c['PD_source_binding']=s['PD_binding']['checked_blocks']==337744 and s['PD_binding']['certificate_sha256']==vc['pd_certificate_sha256'] and s['PD_binding']['checkpoint_sha256']==vc['checkpoint_sha256'] and s['array_sha256']==vc['array_sha256']
# Independent longer rational series (64/80 vs producer32/40).
def atanbox(q):
 lo=sum((F((-1)**m,q**(2*m+1)*(2*m+1)) for m in range(64)),F(0));return lo,lo+F(1,q**129*129)
a,b=atanbox(5);d,e=atanbox(239);piL=16*a-4*e;piU=16*b-4*d
l2L=2*sum((F(1,3**(2*m+1)*(2*m+1)) for m in range(80)),F(0));l2U=l2L+F(2,3**161*161)/(1-F(1,9));pi=s['constant_intervals']['pi'];k1=s['constant_intervals']['k1_log2_over_pi'];c['independent_constants']=Q(pi[0])<=piL<=piU<=Q(pi[1]) and Q(k1[0])<=l2L/piU<=l2U/piL<=Q(k1[1]) and Q(s['constant_intervals']['M_upper'])>=1/piL
with tempfile.TemporaryDirectory(prefix='qcl-endpoint-review-') as td:
 def binary(file,name,start=None,count=None):
  dest=Path(td)/'dataset.bin';args=['h5dump','-d',name]
  if start is not None:args+=['-s',start,'-c',count]
  subprocess.run(args+['-b','LE','-o',str(dest),str(file)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE);return dest.read_bytes()
 def vals(file,name,start=None,count=None):
  buf=binary(file,name,start,count);return struct.unpack('<'+'d'*(len(buf)//8),buf)
 en=o/'endpoint_enclosures.h5';lr=vals(en,'positive_weighted_Gamma_lower_real');li=vals(en,'positive_weighted_Gamma_lower_imag');ur=vals(en,'positive_weighted_Gamma_upper_real');ui=vals(en,'positive_weighted_Gamma_upper_imag');ell=vals(en,'support_left_interval');rr=vals(en,'support_right_interval');A,B=vals(en,'native_window');c['1100_entrywise_complex_boxes']=len(lr)==len(li)==len(ur)==len(ui)==1100 and all(math.isfinite(x) for xs in (lr,li,ur,ui) for x in xs) and all(a<=b for a,b in zip(lr,ur)) and all(a<=b for a,b in zip(li,ui))
 vp=(p/cfg['V4_input']).resolve();base=vp.parent;source=(base/vc['spectral_selected']).resolve();hR=vals(source,'h_H_real');hI=vals(source,'h_H_imag');TR=vals(source,'Tplus_real');TI=vals(source,'Tplus_imag');wk=vals(source,'weight_k')
 interleaved=lambda re,im:b''.join(struct.pack('<dd',a,b) for a,b in zip(re,im));c['small_source_array_hashes']=hashlib.sha256(interleaved(hR,hI)).hexdigest()==vc['array_sha256']['hcell'] and hashlib.sha256(interleaved(TR,TI)).hexdigest()==vc['array_sha256']['Tplus'] and hashlib.sha256(struct.pack('<11d',*wk)).hexdigest()==vc['array_sha256']['weight_momentum']
 gaprows=[]
 for kk in range(11):
  for q in range(4):
   sign=1 if q%2==0 else -1;edge=(ell,rr,(A,A),(B,B))[q];exactrows=[]
   for i in range(45):
    cell,a=divmod(i,5);radius=F(0);diag=None
    for j in range(45):
     cj,bj=divmod(j,5);local=cj==cell
     if local:hidx=(bj*5+a)*11+kk;hr,hi=Q(hR[hidx]),Q(hI[hidx])
     elif cj==cell+1:idx=bj*5+a;hr,hi=Q(TR[idx]),Q(TI[idx])
     elif cj==cell-1:idx=a*5+bj;hr,hi=Q(TR[idx]),-Q(TI[idx])
     else:hr=hi=F(0)
     reL=reU=sign*hr;imL=imU=sign*hi
     if local:
      ix=((bj*5+a)*11+kk)*4+q;reL-=Q(ur[ix]);reU-=Q(lr[ix]);imL-=Q(ui[ix]);imU-=Q(li[ix])
     if i==j:diag=reL-Q(edge[1]) if sign==1 else reL+Q(edge[0])
     else:radius+=max(abs(reL),abs(reU))+max(abs(imL),abs(imU))
    exactrows.append(diag-radius)
   exactmin=min(exactrows);reported=s['per_k'][str(kk+1)]['lower_gaps'][q];assert F(0)<Q(reported)<=exactmin,(kk,q,reported,float(exactmin));gaprows.append({'k':kk+1,'q':q+1,'reported_lower':reported,'exact_box_Gershgorin':float(exactmin),'conservative_slack':float(exactmin-Q(reported))})
 c['all44_exact_rational_gap_reconstructions']=True;r['gap_replay']=gaprows
 # Geometry uses exact dyadic canonical centres generated exactly as V4 Float expressions.
 n=vc['expected_NE'];a=float(vc['PI_origin_dimensionless']);h=float(vc['PI_spacing_dimensionless']);centres=[a+j*h for j in range(n)];qh=Q(h);ellQ=Q(centres[0])-qh;rrQ=Q(centres[-1])+qh
 c['dyadic_support_enclosed']=Q(ell[0])<=ellQ<=Q(ell[1]) and Q(rr[0])<=rrQ<=Q(rr[1]) and ellQ<Q(A)<Q(B)<rrQ;c['adjacent_hats_overlap']=all(0<Q(y)-Q(x)<2*qh for x,y in zip(centres,centres[1:]));c['only_two_endpoint_hats_may_leave_window']=Q(centres[2])-qh>=Q(A) and Q(centres[-3])+qh<=Q(B)
 def area(cen,x):
  if x<=cen-qh:return F(0)
  if x>=cen+qh:return qh
  return (x-cen+qh)**2/(2*qh) if x<cen else qh-(cen+qh-x)**2/(2*qh)
 indices=[0,1,n-2,n-1];areas={j:(area(Q(centres[j]),Q(A)),qh-area(Q(centres[j]),Q(B))) for j in indices}
 gamma=(base/vc['checkpoint']).resolve();endpoint={}
 for start in (0,n-2):
  real=vals(gamma,'Gamma_real',f'0,0,0,{start}','5,5,11,2');imag=vals(gamma,'Gamma_imag',f'0,0,0,{start}','5,5,11,2')
  for shift in (0,1):
   endpoint[start+shift]=[Q(real[((bb*5+aa)*11+k)*2+shift]) for bb in range(5) for aa in range(5) for k in range(11)]
 beta_rows=[]
 for kk in range(11):
  left=right=F(0)
  for j in indices:
   trace=sum((endpoint[j][(aa*5+aa)*11+kk] for aa in range(5)),F(0));assert trace>0
   left+=trace*areas[j][0];right+=trace*areas[j][1]
  g=s['per_k'][str(kk+1)]['lower_gaps'];upper=(left/Q(g[2])**2+right/Q(g[3])**2)/(2*Q(pi[0]));reported=Q(s['per_k'][str(kk+1)]['ghost_mass_operator_upper']);assert reported>=upper>0
  beta_rows.append({'k':kk+1,'exact_dyadic_trace_area_left':float(left),'exact_dyadic_trace_area_right':float(right),'replayed_upper_with_reported_pi_lower':float(upper),'reported':float(reported),'relative_conservative_slack':float((reported-upper)/upper)})
 c['all11_ghost_mass_upper_reconstructions']=True;r['ghost_replay']=beta_rows;r['source_Gamma_hyperslab_scope']={'logical_E_indices_1based':[j+1 for j in indices],'all_k':11,'block_shape':[5,5],'payload_real_plus_imag_bytes':4*11*25*16,'full_Gamma_hash_or_scan_repeated':False}
 sheet=Q(s['units']['spin_degeneracy'])*5*sum((Q(wk[k])*Q(s['per_k'][str(k+1)]['ghost_mass_operator_upper']) for k in range(11)),F(0));c['sheet_upper_exact_reconstructed']=Q(s['ghost_sheet_number_upper'])>=sheet and Q(s['ghost_sheet_over_target_upper'])>=Q(s['ghost_sheet_number_upper'])/Q(s['units']['target_dimensionless'])
 r['geometry']={'support_left':float(ellQ),'support_right':float(rrQ),'native_A':A,'native_B':B,'endpoint_hat_areas':{str(j+1):[float(x),float(y)] for j,(x,y) in areas.items()}}
log=(p/'actual-v2-root.log').read_text();rss=int(re.search(r'Maximum resident set size \(kbytes\): (\d+)',log)[1]);r['cost']={'process_seconds':12.37,'script_seconds':s['elapsed_seconds'],'MaxRSS_KiB':rss,'MaxRSS_bytes':rss*1024,'swap0':'Swaps: 0' in log,'exit0':'Exit status: 0' in log,'output_bytes':sum(f.stat().st_size for f in o.rglob('*') if f.is_file())};c['execution_log']=r['cost']['swap0'] and r['cost']['exit0'] and r['cost']['output_bytes']<16*1024**2;r['minimum_gap']=min(gaprows,key=lambda x:x['reported_lower']);r['all_checks_pass']=all(c.values());assert r['all_checks_pass'],{a:b for a,b in c.items() if not b};Path('.build/campaign133/reviews/endpoint-gap-certificate-result-audit.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'summary_sha256':r['summary_sha256'],'checks':c,'minimum_gap':r['minimum_gap'],'geometry':r['geometry'],'cost':r['cost'],'ghost_sheet_number_upper':s['ghost_sheet_number_upper'],'ghost_sheet_over_target_upper':s['ghost_sheet_over_target_upper']},indent=2))
