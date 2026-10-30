from pathlib import Path
import pandas as pd,json,sys
sys.stdout.reconfigure(encoding='utf-8')
p=Path(__file__).resolve().parent
rows=[];checks=[]
for c in ['HRS','ELSA']:
 d=pd.read_csv(p/f'data_private/{c}_longitudinal.csv',dtype={'exposure':str})
 old=pd.read_csv(p/'audit/attrition_reasons.csv')
 for j in [1,2]:
  for g,x in [('all',d)]+[(g,d[d.exposure.eq(g)]) for g in ['00','01','10','11']]:
   risk=pd.Series(True,index=x.index) if j==1 else x.y1.eq(0)
   missing=risk&x[f'y{j}'].isna()
   dead=missing&x[f'dead{j}'].eq(1) if c=='HRS' else pd.Series(False,index=x.index)
   invalid=missing&~dead&x[f'iwstat{j}'].eq(1)
   other=missing&~dead&~invalid
   assert int(missing.sum())==int(dead.sum()+invalid.sum()+other.sum())
   rows.append(dict(cohort=c,interval=j,group=g,risk_n=int(risk.sum()),missing_total=int(missing.sum()),death_recorded=int(dead.sum()) if c=='HRS' else None,invalid_outcome_time=int(invalid.sum()),other_unknown=int(other.sum())))
  aggregate=next(z for z in rows if z['cohort']==c and z['interval']==j and z['group']=='all')
  src=old[(old.cohort==c)&(old.interval==j)].set_index('reason').n
  for v,reason in [('missing_total','missing_total'),('invalid_outcome_time','interview_without_valid_outcome_or_time'),('other_unknown','other_noninterview_or_unknown')]+([('death_recorded','death_recorded')] if c=='HRS' else []):
   summed=sum(z[v] for z in rows if z['cohort']==c and z['interval']==j and z['group']!='all')
   passed=summed==aggregate[v]==int(src[reason]);assert passed
   checks.append(dict(cohort=c,interval=j,variable=v,group_sum=summed,original_total=int(src[reason]),passed=passed))
pd.DataFrame(rows).to_csv(p/'audit/attrition_by_pattern.csv',index=False)
pd.DataFrame(checks).to_csv(p/'audit/attrition_group_checks.csv',index=False)
print(json.dumps({'rows':len(rows),'checks':len(checks),'all_passed':all(z['passed'] for z in checks)}))