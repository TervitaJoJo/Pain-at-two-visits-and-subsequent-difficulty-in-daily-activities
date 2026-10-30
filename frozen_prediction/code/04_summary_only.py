"""Summarize frozen prediction output; no model refitting or participant export."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
r=root/'results';a=root/'audit'
p=pd.read_csv(r/'performance_all.csv');b=pd.read_csv(r/'two_sample_bootstrap.csv')
keys=['cohort','validation','variant','weighting','replicate','horizon']
x=p[p.model=='M0'].merge(p[p.model=='M1'],on=keys,suffixes=('_M0','_M1'),validate='one_to_one')
for k in ['auc','brier','brier_hajek']:x['delta_'+k]=x[k+'_M1']-x[k+'_M0']
x.to_csv(r/'paired_increment_all.csv',index=False)
bb=b[b.model=='M0'].merge(b[b.model=='M1'],on=['cohort','bootstrap','horizon'],suffixes=('_M0','_M1'),validate='one_to_one')
external=[]
for _,z in x[(x.validation=='frozen_external')&(x.variant=='common')&(x.weighting=='conditional')].iterrows():
 t=bb[(bb.cohort==z.cohort)&(bb.horizon==z.horizon)]
 row=z.to_dict()
 for k in ['auc','brier']:
  d=t[k+'_M1']-t[k+'_M0'];lo,hi=np.quantile(d,[.025,.975]);row['delta_'+k+'_lower']=lo;row['delta_'+k+'_upper']=hi
 row['bootstrap_n']=len(t);external.append(row)
e=pd.DataFrame(external);e.to_csv(r/'external_increment_summary.csv',index=False)
cv=x[x.validation=='local_CV'].groupby(['cohort','variant','weighting','horizon']).agg(auc_M0=('auc_M0','mean'),auc_M1=('auc_M1','mean'),delta_auc=('delta_auc','mean'),delta_auc_min=('delta_auc','min'),delta_auc_max=('delta_auc','max'),brier_M0=('brier_M0','mean'),brier_M1=('brier_M1','mean'),delta_brier=('delta_brier','mean'),repeats=('replicate','count')).reset_index()
cv.to_csv(r/'cv_increment_summary.csv',index=False)

# Aggregate-only check does not require or distribute participant data.
expected=(a/'bootstrap_first_run_sha256.txt').read_text(encoding='utf-8').strip()
assert hashlib.sha256((r/'two_sample_bootstrap.csv').read_bytes()).hexdigest()==expected
assert len(e)==4 and e.bootstrap_n.eq(500).all()
print('Aggregate summaries regenerated; bootstrap hash verified. Participant-level checks were not rerun.')
print(e[e.horizon==2][['cohort','delta_auc','delta_auc_lower','delta_auc_upper']].to_string(index=False))
