import os
DATA_ROOT=__import__("pathlib").Path(os.environ["COHORT_DATA_ROOT"]).resolve()
from pathlib import Path
from datetime import datetime
import json,sys,hashlib
import numpy as np
import pandas as pd
import pyreadstat
sys.stdout.reconfigure(encoding='utf-8')
P=Path(__file__).resolve().parents[1];R=P.parents[1]
checks=[];sources=[];flows=[];summaries=[];missing=[];countrycounts=[];waveaudit=[]
ADL=['dressa','batha','eata','beda','toilta'];COV=['age','female','education','smoking','arthritis','cvd','cesd']
def check(c,n,ok,detail=''):
 checks.append(dict(cohort=c,check=n,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError((c,n,detail))
def binary(v):return v.where(v.isin([0,1]))
def anybin(d):return pd.Series(np.where(d.eq(1).any(axis=1),1,np.where(d.eq(0).all(axis=1),0,np.nan)),index=d.index)
def read(p,cols,idcol,c):
 d,m=pyreadstat.read_dta(p,usecols=list(dict.fromkeys([idcol]+cols)))
 check(c,p.name+'_unique_id',d[idcol].notna().all() and not d[idcol].duplicated().any())
 sources.append(dict(cohort=c,path=str(p),bytes=p.stat().st_size,mtime=p.stat().st_mtime,columns=list(d.columns)))
 d.index=d[idcol].astype(str);return d
# CHARLS harmonized variables and raw pain/response provenance.
c='CHARLS';p=DATA_ROOT/'CHARLS/Harmonized_CHARLS/H_CHARLS_D_Data.dta'
sufs=ADL+['iwstat','agey','iwy','iwm','diabe','arthre','hearte','stroke','smoken','cesd10','cesd10m']
d=read(p,['ID_w1','ragender','raeducl']+[f'r{w}{s}' for w in range(1,5) for s in sufs],'ID',c)
W={}
for w in range(1,5):
 z=pd.DataFrame(index=d.index);v=lambda s:d[f'r{w}{s}']
 z['age']=v('agey').where(v('agey').between(0,120));z['female']=d.ragender.map({1:0.,2:1.})
 z['education']=d.raeducl.map({0:1,1:1,2:2,3:3});z['smoking']=binary(v('smoken'));z['arthritis']=binary(v('arthre'))
 z['cvd']=anybin(d[[f'r{w}hearte',f'r{w}stroke']]);z['cesd']=v('cesd10').where(v('cesd10m').eq(0)&v('cesd10').between(0,30))
 aa=d[[f'r{w}{s}' for s in ADL]];check(c,f'w{w}_ADL_codes',aa.stack().dropna().isin([0,1]).all());z['adl']=anybin(aa)
 z['month']=(v('iwy')*12+v('iwm')).where(v('iwm').between(1,12)&v('iwy').between(2010,2022))
 z['interview']=v('iwstat').eq(1);z['diabetes']=binary(v('diabe'));W[w]=z
for w,y in [(1,2011),(2,2013)]:
 p=next((DATA_ROOT/'CHARLS'/str(y)).glob('*[Hh]ealth*[Ff]unctioning*.dta'))
 raw=read(p,['da041','proxy'] if w==1 else ['wb16','db032'],'ID',c)
 key=d.ID_w1.astype(str) if w==1 else pd.Series(d.index,index=d.index)
 check(c,f'w{w}_ID_coverage',key[d[f'r{w}iwstat'].eq(1)].isin(raw.index).mean()>.99)
 def rv(n):return key.map(raw[n])
 W[w]['pain']=rv('da041').map({1:1.,2:0.}) if w==1 else rv('wb16').map({1:0.,2:1.,3:1.,4:1.,5:1.})
 W[w]['self']=rv('proxy').eq(0) if w==1 else rv('db032').isin([1,2,3])
 W[w]['self_only']=rv('proxy').eq(0) if w==1 else rv('db032').eq(1)
COHORTS={'CHARLS':(W,[(1,2,3,4)])}
# SHARE direct raw release fields; never use other studies' constructed working datasets.
c='SHARE';W={}
for w in range(5,10):
 def mod(su,cols):
  p=next((DATA_ROOT/'SHARE').rglob(f'sharew{w}_rel9-0-0_{su}.dta'));return read(p,cols,'mergeid',c)
 phcols=['country','ph084_','ph054_','ph006d1','ph006d4','ph006d5','ph006d19','ph006d20','ph048d3','ph048d5']+[f'ph049d{i}' for i in range(1,7)]
 ph=mod('ph',phcols);cv=mod('cv_r',['age_int','gender','int_year','int_month']).reindex(ph.index)
 z=pd.DataFrame(index=ph.index);z['country']=ph.country;z['age']=cv.age_int.where(cv.age_int.between(0,120));z['female']=cv.gender.map({1:0.,2:1.})
 z['month']=(12*cv.int_year+cv.int_month).where(cv.int_year.between(2000,2026)&cv.int_month.between(1,12));z['interview']=z.month.notna()
 z['adl']=anybin(ph[[f'ph049d{i}' for i in [1,3,4,5,6]]].map(lambda x:x if x in [0,1] else np.nan))
 z['function3']=ph[['ph049d2','ph048d3','ph048d5']].map(lambda x:x if x in [0,1] else np.nan).sum(axis=1,min_count=3)
 z['pain']=ph.ph084_.map({1:1.,5:0.});z['self']=ph.ph054_.isin([1,2]);z['self_only']=ph.ph054_.eq(1)
 z['diabetes']=binary(ph.ph006d5);z['arthritis']=anybin(ph[['ph006d19','ph006d20']].map(lambda x:x if x in [0,1] else np.nan));z['cvd']=anybin(ph[['ph006d1','ph006d4']].map(lambda x:x if x in [0,1] else np.nan))
 if w in [5,6]:
  ed=mod('gv_isced',['isced1997_r']).reindex(ph.index);br=mod('br',['br001_','br002_']+(['br022_'] if w==5 else [])).reindex(ph.index);mh=mod('gv_health',['eurod']+[f'euro{i}' for i in range(1,13)]).reindex(ph.index)
  z['education']=ed.isced1997_r.map({0:1,1:1,2:1,3:2,4:2,5:3,6:3})
  z['smoking']=br.br002_.map({1:1.,5:0.});z.loc[br.br001_.eq(5),'smoking']=0.
  if w==5:
   z.loc[br.br022_.isin([1,2]),'smoking']=0.;z.loc[br.br022_.eq(3),'smoking']=1.
  allcomp=mh[[f'euro{i}' for i in range(1,13)]].isin([0,1]).all(axis=1)
  z['cesd']=mh.eurod.where(mh.eurod.between(0,12)&allcomp)
 W[w]=z
 waveaudit.append(dict(cohort=c,wave=w,n=len(z),known_pain=int(z.pain.notna().sum()),known_adl=int(z.adl.notna().sum()),respondent_involved=int(z.self.sum()),date_min=float(z.month.min()),date_max=float(z.month.max())))
idx=W[5].index
for w in range(6,10):idx=idx.union(W[w].index)
W={w:z.reindex(idx) for w,z in W.items()};COHORTS['SHARE']=(W,[(5,6,7,8),(6,7,8,9)])
for c,(W,windows) in COHORTS.items():
 assigned=pd.Series(False,index=next(iter(W.values())).index);parts=[]
 for a,b,k,l in windows:
  A,B=W[a],W[b]
  def stage(n,m):flows.append(dict(cohort=c,window=a,stage=n,n=int(m.sum())))
  mask=A.interview.eq(True)&A.age.ge(50);stage('age_and_interview',mask)
  mask &= A.diabetes.eq(1);stage('diabetes',mask)
  mask &= A.self.eq(True)&B.self.eq(True)&B.interview.eq(True);stage('two_self_interviews',mask)
  delta=B.month-A.month;mask &= delta.between(12,36);stage('12_36_month_exposure_interval',mask)
  mask &= A.adl.eq(0)&B.adl.eq(0);stage('both_ADL_free',mask)
  mask &= A.pain.notna()&B.pain.notna();stage('both_pain_known',mask)
  take=mask&~assigned;assigned|=take;stage('earliest_eligible',take)
  keep=COV+(['country','function3'] if c=='SHARE' else [])
  z=A.loc[take,keep].copy();z['entry']=a;z['exposure']=A.loc[take,'pain'].astype(int).astype(str)+B.loc[take,'pain'].astype(int).astype(str)
  z['exposure_years']=delta.loc[take]/12;z['landmark_month']=B.loc[take,'month'];z['self_only']= (A.self_only.eq(True)&B.self_only.eq(True)).loc[take].astype(int)
  for j,w in enumerate([k,l],1):
   F=W[w];valid=F.interview.eq(True)&F.month.gt(B.month)
   z[f'y{j}']=F.adl.where(valid).loc[take];z[f't{j}']=((F.month-B.month)/12).where(valid).loc[take]
  parts.append(z)
 z=pd.concat(parts);check(c,'unique_eligible_ids',not z.index.duplicated().any());check(c,'known_exposure',z.exposure.isin(['00','01','10','11']).all())
 check(c,'positive_t1',z.loc[z.y1.notna(),'t1'].gt(0).all());check(c,'ordered_t2',z.loc[z.y1.eq(0)&z.y2.notna(),'t2'].gt(z.loc[z.y1.eq(0)&z.y2.notna(),'t1']).all())
 z.insert(0,'pid',range(1,len(z)+1));z.reset_index(drop=True,inplace=True)
 z.to_csv(P/f'data_private/{c}_longitudinal.csv',index=False)
 if c=='SHARE':
  z.loc[z.entry.eq(5)].to_csv(P/'data_private/SHARE_primary_longitudinal.csv',index=False)
 for var in COV+(['function3'] if c=='SHARE' else []):missing.append(dict(cohort=c,variable=var,n=len(z),missing=int(z[var].isna().sum())))
 for group,s in z.groupby('exposure'):
  obs=s.y1.notna();ev=s.y1.eq(1)|(s.y1.eq(0)&s.y2.eq(1));use2=s.y1.eq(0)&s.y2.notna()
  summaries.append(dict(cohort=c,group=group,eligible=len(s),with_followup=int(obs.sum()),events=int(ev.sum()),intervals=int(obs.sum()+use2.sum()),observed_years=s.loc[obs,'t1'].sum()+(s.loc[use2,'t2']-s.loc[use2,'t1']).sum()))
 if c=='SHARE':
  for country,s in z.groupby('country'):countrycounts.append(dict(country=country,eligible=len(s),observed=int(s.y1.notna().sum()),events=int((s.y1.eq(1)|(s.y1.eq(0)&s.y2.eq(1))).sum())))
 print(c,'eligible',len(z),'followup',z.y1.notna().sum(),'events',(z.y1.eq(1)|(z.y1.eq(0)&z.y2.eq(1))).sum(),'entries',z.entry.value_counts().to_dict(),flush=True)
for name,rows in [('extraction_checks',checks),('sample_flow',flows),('group_counts',summaries),('covariate_missingness',missing),('SHARE_countries',countrycounts),('SHARE_wave_availability',waveaudit)]:pd.DataFrame(rows).to_csv(P/f'audit/{name}.csv',index=False)
(P/'audit/data_sources.json').write_text(json.dumps(sources,indent=2),encoding='utf-8')
(P/'audit/extraction_manifest.json').write_text(json.dumps({'timestamp':datetime.now().astimezone().isoformat(),'all_checks_passed':all(x['passed'] for x in checks),'checks':len(checks),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),encoding='utf-8')



