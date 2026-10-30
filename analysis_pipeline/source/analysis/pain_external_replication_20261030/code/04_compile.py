from pathlib import Path
from datetime import datetime
import json,hashlib
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[1];R=P.parents[1]
old=R/'analysis/pain_transition_diabetes_adl_20261029'
rows=[];group_rows=[];checks=[]
for c in ['HRS','ELSA','SHARE','CHARLS']:
 if c in ['HRS','ELSA']:
  est=pd.read_csv(old/f'results/{c}_contrasts.csv');src=old/f'data_private/{c}_longitudinal.csv'
 else:
  est=pd.read_csv(P/f'results/{c}_estimates.csv');src=P/'data_private'/f'{c}_{"primary_" if c=="SHARE" else ""}longitudinal.csv'
 d=pd.read_csv(src,dtype={'exposure':str});events=d.y1.eq(1)|(d.y1.eq(0)&d.y2.eq(1));obs=d.y1.notna()
 e=est[(est.variant=='primary_MI')&(est.contrast=='11_vs_00')].iloc[0]
 assert int(e.persons)==int(obs.sum())
 rows.append(dict(cohort=c,eligible=len(d),with_followup=int(obs.sum()),events=int(events.sum()),aHR=float(e.estimate),lower=float(e.lower),upper=float(e.upper),p=float(e.p),classification='Measurement-adapted exploratory' if c=='CHARLS' else 'Original cohort' if c in ['HRS','ELSA'] else 'Additional conceptual replication'))
 for g,s in d.groupby('exposure'):
  ev=s.y1.eq(1)|(s.y1.eq(0)&s.y2.eq(1));o=s.y1.notna();n2=s.y1.eq(0)&s.y2.notna()
  group_rows.append(dict(cohort=c,group=g,eligible=len(s),with_followup=int(o.sum()),events=int(ev.sum()),intervals=int(o.sum()+n2.sum()),observed_years=float(s.loc[o,'t1'].sum()+(s.loc[n2,'t2']-s.loc[n2,'t1']).sum())))
 if c in ['SHARE','CHARLS']:
  diag=pd.read_csv(P/f'audit/{c}_model_diagnostics.csv');main=diag[diag.variant=='primary_MI']
  assert len(main)==50 and main.converged.all() and main.warnings.sum()==0
  assert (main.persons==obs.sum()).all() and (main.events==events.sum()).all()
  checks.append({'cohort':c,'main_models_converged':50,'all_model_fits':len(diag),'model_warnings':int(diag.warnings.sum()),'counts_reconciled':True})
pd.DataFrame(rows).to_csv(P/'results/four_cohort_primary_comparison.csv',index=False)
pd.DataFrame(group_rows).to_csv(P/'results/four_cohort_group_counts.csv',index=False)
v=json.loads((P/'audit/independent_model_verification.json').read_text())
assert len(v['checks'])==6 and all(x['passed'] for x in v['checks'])
files=[]
for p in P.rglob('*'):
 if p.is_file() and p.suffix in ['.py','.R','.md','.csv','.json'] and p.name!='analysis_manifest.json':
  files.append({'path':str(p.relative_to(P)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(P/'audit/analysis_manifest.json').write_text(json.dumps({'timestamp':datetime.now().astimezone().isoformat(),'checks':checks,'independent_model_checks':6,'files':files},indent=2),encoding='utf-8')
print(pd.DataFrame(rows).to_string(index=False))

