from pathlib import Path
import pandas as pd,json,sys
sys.stdout.reconfigure(encoding='utf-8')
p=Path(__file__).parent
rows=[];cells=[];proxies=[];attr=[]
for c in ['HRS','ELSA']:
 d=pd.read_csv(p/f'data_private/{c}_longitudinal.csv',dtype={'exposure':str})
 for scope,x in [('all',d)]+[(str(w),z) for w,z in d.groupby('entry')]:
  for v in ['diabetes_treatment','landmark_function3']+(['bmi','diabetes_duration'] if c=='HRS' else []):
   q=x[v].dropna(); rows.append(dict(cohort=c,window=scope,variable=v,n=len(x),missing=x[v].isna().sum(),min=q.min(),max=q.max(),mean=q.mean(),sd=q.std()))
 for g,x in d.groupby('exposure'):
  for v in ['diabetes_treatment','landmark_function3']:
   for lev in sorted(x[v].dropna().unique()):cells.append(dict(cohort=c,group=g,variable=v,level=lev,n=int(x[v].eq(lev).sum())))
  pids=set(); ni=ne=0
  for j in [1,2]:
   risk=pd.Series(True,index=x.index) if j==1 else x.y1.eq(0)
   valid=risk&x[f'y{j}'].notna();prox=valid&x[f'proxy{j}'].eq(1)
   pids.update(x.loc[prox,'pid']);ni+=prox.sum();ne+=x.loc[prox,f'y{j}'].sum()
  proxies.append(dict(cohort=c,group=g,eligible=len(x),proxy_persons=len(pids),proxy_intervals=int(ni),proxy_events=int(ne)))
 for j in [1,2]:
  risk=pd.Series(True,index=d.index) if j==1 else d.y1.eq(0)
  miss=risk&d[f'y{j}'].isna()
  dead=miss&d[f'dead{j}'].eq(1) if c=='HRS' else pd.Series(False,index=d.index)
  interv=miss&~dead&d[f'iwstat{j}'].eq(1)
  other=miss&~dead&~interv
  assert int(miss.sum())==int(dead.sum()+interv.sum()+other.sum())
  for reason,mask in [('missing_total',miss),('death_recorded' if c=='HRS' else 'death_not_separately_available',dead),('interview_without_valid_outcome_or_time',interv),('other_noninterview_or_unknown',other)]:attr.append(dict(cohort=c,interval=j,reason=reason,n=int(mask.sum())))
for name,records in [('new_variable_audit',rows),('new_variable_cells',cells),('proxy_by_exposure',proxies),('attrition_reasons',attr)]:pd.DataFrame(records).to_csv(p/f'audit/{name}.csv',index=False)
print(pd.DataFrame(rows).query("window=='all'").to_string(index=False))
print(pd.DataFrame(proxies).groupby('cohort')[['proxy_persons','proxy_intervals','proxy_events']].sum().to_string())
