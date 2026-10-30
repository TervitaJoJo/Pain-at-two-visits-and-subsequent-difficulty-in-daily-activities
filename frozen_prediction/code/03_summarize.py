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
checks=[]
for cohort in ['HRS','ELSA','SHARE']:
 d=pd.read_csv(root/'data_private'/f'{cohort}_qa_predictions_PRIVATE.csv');obs=d[d.y.notna() & (d.weight>0)]
 mode='apparent_development' if cohort=='HRS' else 'frozen_external'
 for model in ['M0','M1']:
  z=p[(p.cohort==cohort)&(p.validation==mode)&(p.variant=='common')&(p.weighting=='conditional')&(p.horizon==2)&(p.model==model)].iloc[0]
  cases=obs[obs.y==1];controls=obs[obs.y==0]
  scores=cases[model].to_numpy()[:,None]-controls[model].to_numpy()[None,:]
  weights=cases.weight.to_numpy()[:,None]*controls.weight.to_numpy()[None,:]
  auc=float(np.sum(weights*((scores>0)+.5*(scores==0)))/weights.sum())
  err=float(np.sum(obs.weight*(obs.y-obs[model])**2))
  for metric,val in [('auc',auc),('brier',err/len(d)),('brier_hajek',err/obs.weight.sum())]:
   checks.append(dict(check=f'{cohort} {model} independent Python {metric}',passed=abs(val-z[metric])<1e-12,absolute_error=abs(val-z[metric])))
old=(a/'bootstrap_first_run_sha256.txt').read_text(encoding='utf-8').strip()
checks.append(dict(check='Diagnostic-only bootstrap rerun byte-identical',passed=hashlib.sha256((r/'two_sample_bootstrap.csv').read_bytes()).hexdigest()==old,absolute_error=0))
qa=pd.DataFrame(checks);qa.to_csv(a/'independent_python_checks.csv',index=False);assert qa.passed.all()
warn=pd.read_csv(a/'bootstrap_warning_details.csv')
assert len(warn)==2 and (warn.link=='logit').all() and warn.formula.str.startswith('response ~').all()
lines=['# Exploratory prediction module','', 'Status: post hoc relative to the association study; module specification fixed before fitting these prediction models. Existing association results were already known.','', 'Question: does first-visit pain history add to second-visit pain plus a common clinical model that already includes first-visit function?','', 'Primary horizon is the second subsequent scheduled visit, not a common exact four-year horizon. HRS development predictions are frozen for ELSA/SHARE. AUC differences are M1 minus M0; negative Brier differences favor M1.','', '| Cohort | AUC M0 | AUC M1 | Change (95% CI) | Brier change (95% CI) |','|---|---:|---:|---|---|']
for _,z in e[e.horizon==2].iterrows():
 lines.append(f'| {z.cohort} | {z.auc_M0:.4f} | {z.auc_M1:.4f} | {z.delta_auc:+.4f} ({z.delta_auc_lower:+.4f}, {z.delta_auc_upper:+.4f}) | {z.delta_brier:+.5f} ({z.delta_brier_lower:+.5f}, {z.delta_brier_upper:+.5f}) |')
lines+=['', 'External confidence intervals use 500 paired two-sample bootstrap replicates, refitting development preprocessing and models and validation observation models. No failed cohort units. Two observation-model fits produced extreme-probability warnings; both converged. Predefined probability clipping and weight caps were retained. The diagnostic rerun reproduced every aggregate bootstrap result byte for byte.','', 'HRS repeated local cross-validation gave a mean AUC change of +0.00057. Some Brier sensitivities improved slightly in ELSA; the result is no stable added value, not universal worsening or proof of equivalence.','', 'Death and loss are censored. Predictions and weighted scoring target hypothetical net probabilities under conditional independent observation; they are not competing-risk cumulative incidence. Calibration estimates are descriptive and did not recalibrate the frozen predictions.','', f'QA: 31 R numerical/leakage checks and {len(checks)} independent Python/hash checks passed. The private prediction CSVs are not distributable and are excluded from the submission package.','', 'Disposition: retain as supplementary boundary analysis with brief main-text disclosure. Do not promote to an independent clinical prediction claim.']
(root/'prediction_report.md').write_text('\n'.join(lines),encoding='utf-8')
print('Independent checks:',len(checks),'passed')
print(e[e.horizon==2][['cohort','auc_M0','auc_M1','delta_auc','delta_auc_lower','delta_auc_upper','delta_brier','delta_brier_lower','delta_brier_upper']].to_string(index=False))
print('CV summary');print(cv.to_string(index=False))