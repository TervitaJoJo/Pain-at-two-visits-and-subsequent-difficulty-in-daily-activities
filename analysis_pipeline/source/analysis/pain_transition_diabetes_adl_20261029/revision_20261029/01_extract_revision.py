import os
DATA_ROOT=__import__("pathlib").Path(os.environ["COHORT_DATA_ROOT"]).resolve()
from pathlib import Path
from datetime import datetime
import json,re,zipfile,hashlib,sys
import numpy as np
import pandas as pd
import pyreadstat
sys.stdout.reconfigure(encoding='utf-8')
ORIGINAL=Path(__file__).resolve().parent.parent; P=ORIGINAL/'revision_20261029'; ROOT=ORIGINAL.parents[1]
CM=json.loads((ROOT/'analysis/second_line_rescreen_20261028/local/community_metadata.json').read_text(encoding='utf-8'))
HM=json.loads((ROOT/'analysis/cross_system_multimorbidity_mortality_20261028/results/metadata_HRS_rand.json').read_text(encoding='utf-8'))
WINS={'HRS':[[10,11,12,13],[11,12,13,14]],'ELSA':[[5,6,7,8],[6,7,8,9]]}
ADL=['dressa','batha','eata','beda','toilta']; FUNC=['walkra','chaira','clim1a']
COV=['age','female','education','smoking','arthritis','cvd','cesd']
checks=[];sources=[];flows=[];summ=[];miss=[];deathaudit=[];timings=[];cells=[];baselines=[]
def check(cohort,name,ok,detail=''):
 checks.append(dict(cohort=cohort,check=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(f'{cohort}: {name}: {detail}')
def binary(v):return v.where(v.isin([0,1]))
def any_binary(a):return pd.Series(np.where(a.eq(1).any(axis=1),1,np.where(a.eq(0).all(axis=1),0,np.nan)),index=a.index)
def core(year):
 stem=f'h{str(year)[2:]}';folder=DATA_ROOT/f'HRS/public_survey_data/{year}'
 dic=folder/f'{stem}sta/{stem.upper()}C_R.dct'
 specs={name:(int(start)-1,int(start)-1+int(width)) for start,typ,name,width in re.findall(r'_column\((\d+)\)\s+(\w+)\s+(\w+)\s+%(\d+)(?:\.\d+)?[sf]',dic.read_text(encoding='utf-8',errors='replace')) if name in ['HHID','PN'] or re.match(r'^[A-Z]C(104|105)$',name)}
 check('HRS',f'{year}_core_fields',len(specs)==4)
 archive=folder/f'{stem}da.zip'
 with zipfile.ZipFile(archive) as z:
  member=next(n for n in z.namelist() if n.upper()==f'{stem.upper()}C_R.DA')
  with z.open(member) as f:raw=pd.read_fwf(f,colspecs=list(specs.values()),names=list(specs),dtype=str)
 key=pd.to_numeric(raw.HHID,errors='coerce')*1000+pd.to_numeric(raw.PN,errors='coerce')
 check('HRS',f'{year}_unique_core_id',key.notna().all() and not key.duplicated().any())
 raw.index=key.astype('int64')
 pa=pd.to_numeric(raw[next(n for n in specs if n.endswith('C104'))],errors='coerce')
 se=pd.to_numeric(raw[next(n for n in specs if n.endswith('C105'))],errors='coerce')
 check('HRS',f'{year}_pain_codes',pa.dropna().isin([-8,1,5,8,9]).all())
 check('HRS',f'{year}_severity_codes',se.dropna().isin([-8,1,2,3,8,9]).all())
 sources.append(dict(path=str(archive),bytes=archive.stat().st_size,member=member,dictionary=str(dic)))
 return pd.DataFrame({'pain':pa.map({1:1.,5:0.}),'severity':se.where(se.isin([1,2,3])).where(pa.eq(1),0).where(pa.isin([1,5]))})
for cohort in ['HRS','ELSA']:
 windows=WINS[cohort];waves=sorted(set(sum(windows,[])))
 ishrs=cohort=='HRS'; allcols={x['name']:x['label'] for x in HM['columns']} if ishrs else CM['ELSA_harmonized']['labels']
 path=HM['path'] if ishrs else CM['ELSA_harmonized']['path']; idcol='hhidpn' if ishrs else 'idauniq'
 suff=ADL+FUNC+['iwstat','proxy','diabe','cesd','cesdm','arthre','hearte','stroke','smoken']+(['agey_m','iwendy','iwendm'] if ishrs else ['agey','iwindy','iwindm','painfr','painlv'])
 cols=[idcol,'ragender','raeduc' if ishrs else 'raeducl']+[f'r{w}{s}' for w in waves for s in suff]
 extra_suff=['bmi','dborlmed','dbinsuln','dbdiagyr'] if ishrs else ['rxdiabi','rxdiabo','rxdiab']
 cols += [f'r{w}{s}' for w in sorted(set(x for win in windows for x in win[:2])) for s in extra_suff]
 if ishrs:cols+=['radyear','radmonth']
 check(cohort,'all_requested_columns_exist',set(cols)<=set(allcols))
 d,metadata=pyreadstat.read_dta(path,usecols=cols,apply_value_formats=False,formats_as_category=False)
 check(cohort,'unique_harmonized_ids',d[idcol].notna().all() and not d[idcol].duplicated().any())
 d.index=d[idcol]
 sources.append(dict(cohort=cohort,path=path,bytes=Path(path).stat().st_size,mtime=Path(path).stat().st_mtime,columns=cols))
 def v(w,s):return pd.to_numeric(d[f'r{w}{s}'],errors='coerce')
 def month(w):
  m=v(w,'iwendm' if ishrs else 'iwindm');y=v(w,'iwendy' if ishrs else 'iwindy')
  return (12*y+m).where(m.between(1,12)&y.between(1990,2024))
 def age(w):return v(w,'agey_m' if ishrs else 'agey')
 adls={}
 for w in waves:
  a=d[[f'r{w}{s}' for s in ADL]]
  check(cohort,f'w{w}_binary_ADL',a.stack().dropna().isin([0,1]).all())
  adls[w]=any_binary(a)
  check(cohort,f'w{w}_no_status2',not v(w,'iwstat').eq(2).any(),int(v(w,'iwstat').eq(2).sum()))
 if ishrs:
  raw={w:core(1990+2*w) for w in [10,11,12]}
  pain={w:raw[w].pain.reindex(d.index) for w in raw};sev={w:raw[w].severity.reindex(d.index) for w in raw}
 else:
  pain={w:binary(v(w,'painfr')) for w in waves};sev={w:v(w,'painlv').where(v(w,'painlv').isin([0,1,2,3])) for w in waves}
 def baseline(w):
  z=pd.DataFrame(index=d.index);z['age']=age(w)
  z['female']=d.ragender.map({1:0.,2:1.})
  z['education']=d['raeduc'].map({1:1,2:2,3:2,4:2,5:3}) if ishrs else d['raeducl'].where(d['raeducl'].isin([1,2,3]))
  z['smoking']=binary(v(w,'smoken'));z['arthritis']=binary(v(w,'arthre'))
  z['cvd']=any_binary(d[[f'r{w}hearte',f'r{w}stroke']])
  z['cesd']=v(w,'cesd').where(v(w,'cesdm').eq(0)&v(w,'cesd').between(0,8))
  f=d[[f'r{w}{s}' for s in FUNC]]
  check(cohort,f'w{w}_binary_function3',f.stack().dropna().isin([0,1]).all())
  z['function3']=f.sum(axis=1,min_count=3);z['pain_severity']=sev[w]
  z['entry']=2010 if w==windows[0][0] else 2012
  oral=binary(v(w,'dborlmed' if ishrs else 'rxdiabo'));ins=binary(v(w,'dbinsuln' if ishrs else 'rxdiabi'))
  z['diabetes_treatment']=np.select([ins.eq(1),ins.eq(0)&oral.eq(1),ins.eq(0)&oral.eq(0)],[2.,1.,0.],default=np.nan)
  if ishrs:
   z['bmi']=v(w,'bmi').where(v(w,'bmi').gt(0))
   duration=v(w,'iwendy')-v(w,'dbdiagyr')
   z['diabetes_duration']=duration.where(duration.ge(0)&duration.le(age(w)+1))
  return z
 assigned=pd.Series(False,index=d.index);parts=[]
 for a,b,c,f in windows:
  def step(name,mask):flows.append(dict(cohort=cohort,sample='longitudinal',window=a,stage=name,n=int(mask.sum())))
  mask=v(a,'iwstat').eq(1)&age(a).ge(50);step('age_interview',mask)
  mask&=v(a,'diabe').eq(1);step('diabetes',mask)
  mask&=v(a,'proxy').eq(0)&v(b,'proxy').eq(0)&v(b,'iwstat').eq(1);step('two_self_interviews',mask)
  mask&=(month(b)-month(a)).between(12,36);step('12_36_month_interval',mask)
  mask&=adls[a].eq(0)&adls[b].eq(0);step('both_ADL_free',mask)
  mask&=pain[a].notna()&pain[b].notna();step('both_pain_known',mask)
  take=mask&~assigned;assigned|=take;step('earliest_eligible',take)
  z=baseline(a).loc[take].copy();z['exposure']=pain[a].loc[take].astype(int).astype(str)+pain[b].loc[take].astype(int).astype(str)
  z['landmark_function3']=d[[f'r{b}{s}' for s in FUNC]].sum(axis=1,min_count=3).loc[take]
  z['landmark_month']=month(b).loc[take];z['exposure_years']=((month(b)-month(a))/12).loc[take]
  for j,w in enumerate([c,f],1):
   obs=v(w,'iwstat').eq(1)&month(w).gt(month(b))
   z[f'y{j}']=adls[w].where(obs).loc[take]
   z[f't{j}']=((month(w)-month(b))/12).where(obs).loc[take]
   z[f'proxy{j}']=v(w,'proxy').loc[take]
   z[f'iwstat{j}']=v(w,'iwstat').loc[take]
   z[f'dead{j}']=v(w,'iwstat').isin([3,5,6]).loc[take].astype(int) if ishrs else np.nan
  if ishrs:
   z['death_month']=(12*d.radyear+d.radmonth).where(d.radmonth.between(1,12)).loc[take]
  parts.append(z)
 longitudinal=pd.concat(parts)
 check(cohort,'unique_longitudinal',not longitudinal.index.duplicated().any())
 if ishrs:
  known=longitudinal[['dead1','dead2']].eq(1).any(axis=1)
  dd=longitudinal.death_month-longitudinal.landmark_month
  check(cohort,'all_known_death_dates_after_landmark',not (known&dd.le(0)).any())
  deathaudit.append(dict(cohort=cohort,known_deaths=int(known.sum()),death_month_missing=int((known&dd.isna()).sum()),min_months=float(dd[known].min()),max_months=float(dd[known].max())))
  # No status-based dead participant may also supply a living ADL report in that wave.
  for j in [1,2]:check(cohort,f'dead_no_live_ADL_{j}',not (longitudinal[f'dead{j}'].eq(1)&longitudinal[f'y{j}'].notna()).any())
 check(cohort,'followup_time_order',not (longitudinal.y1.notna()&longitudinal.y2.notna()&longitudinal.t2.le(longitudinal.t1)).any())
 # Cross-sectional entry does not depend on any subsequent measurement.
 assigned=pd.Series(False,index=d.index);parts=[]
 for a,*_ in windows:
  mask=v(a,'iwstat').eq(1)&v(a,'proxy').eq(0)&age(a).ge(50)&v(a,'diabe').eq(1)&pain[a].notna()&adls[a].notna()
  take=mask&~assigned;assigned|=take
  flows.append(dict(cohort=cohort,sample='cross_sectional',window=a,stage='earliest_eligible',n=int(take.sum())))
  z=baseline(a).loc[take].copy();z['pain']=pain[a].loc[take];z['adl']=adls[a].loc[take];parts.append(z)
 cross=pd.concat(parts);check(cohort,'unique_cross_sectional',not cross.index.duplicated().any())
 for sample,z in [('longitudinal',longitudinal),('cross_sectional',cross)]:
  z=z.copy();z.insert(0,'pid',range(1,len(z)+1));z.reset_index(drop=True,inplace=True)
  original=pd.read_csv(ORIGINAL/f'data_private/{cohort}_{sample}.csv',dtype={'exposure':str})
  compare=z[original.columns].reset_index(drop=True).copy()
  for col in original.columns:
   if col!='exposure':compare[col]=pd.to_numeric(compare[col],errors='raise').astype(float)
  pd.testing.assert_frame_equal(compare,original,check_dtype=False,check_exact=False,atol=1e-9,rtol=1e-10)
  check(cohort,f'{sample}_all_original_columns_identical',True)
  z.to_csv(P/f'data_private/{cohort}_{sample}.csv',index=False,encoding='utf-8')
  for col in COV+['function3','pain_severity']:
   miss.append(dict(cohort=cohort,sample=sample,variable=col,n=len(z),missing=int(z[col].isna().sum())))
  miss.append(dict(cohort=cohort,sample=sample,variable='any_main_covariate',n=len(z),missing=int(z[COV].isna().any(axis=1).sum())))
  groupcol='exposure' if sample=='longitudinal' else 'pain'
  for g,sub in z.groupby(groupcol):
   for variable in COV+['function3','pain_severity']:
    x=sub[variable].dropna()
    baselines.append(dict(cohort=cohort,sample=sample,group=g,variable=variable,n=len(sub),known=len(x),mean=x.mean(),sd=x.std(),median=x.median(),q25=x.quantile(.25),q75=x.quantile(.75)))
    if variable in ['female','education','smoking','arthritis','cvd','pain_severity']:
     for level,count in sub[variable].value_counts(dropna=False).items():cells.append(dict(cohort=cohort,sample=sample,group=g,variable=variable,level=level,count=int(count)))
   if sample=='longitudinal':
    observed=sub.y1.notna();event=sub.y1.eq(1)|(sub.y1.eq(0)&sub.y2.eq(1))
    keep2=sub.y1.eq(0)&sub.y2.notna()
    years=sub.loc[observed,'t1'].sum()+(sub.loc[keep2,'t2']-sub.loc[keep2,'t1']).sum()
    summ.append(dict(cohort=cohort,sample=sample,group=g,eligible=len(sub),with_followup=int(observed.sum()),events=int(event.sum()),intervals=int(observed.sum()+keep2.sum()),observed_interval_years=years,no_first_followup=int((~observed).sum()),missing_second_after_negative_first=int((sub.y1.eq(0)&sub.y2.isna()).sum())))
   else:summ.append(dict(cohort=cohort,sample=sample,group=g,eligible=len(sub),with_followup=len(sub),events=int(sub.adl.sum())))
  if sample=='longitudinal':
   for col in ['exposure_years','t1','t2']:
    x=z[col].dropna();timings.append(dict(cohort=cohort,variable=col,n=len(x),min=x.min(),q25=x.quantile(.25),median=x.median(),q75=x.quantile(.75),max=x.max()))
 expected= {'HRS':dict(n=2993,events=510,first=2650,groups={'00':1670,'01':349,'10':307,'11':667}), 'ELSA':dict(n=643,events=99,first=544,groups={'00':298,'01':94,'10':81,'11':170})}[cohort]
 check(cohort,'screening_n_reproduced',len(longitudinal)==expected['n'])
 check(cohort,'screening_groups_reproduced',longitudinal.exposure.value_counts().to_dict()==expected['groups'])
 check(cohort,'screening_first_observed_reproduced',longitudinal.y1.notna().sum()==expected['first'])
 check(cohort,'screening_events_reproduced',(longitudinal.y1.eq(1)|(longitudinal.y1.eq(0)&longitudinal.y2.eq(1))).sum()==expected['events'])
 print(cohort,'landmark',len(longitudinal),'cross_sectional',len(cross),'capacity reconciliation passed',flush=True)
for name,rows in [('extraction_checks',checks),('sample_flow',flows),('group_outcomes',summ),('missingness',miss),('death_audit',deathaudit),('timing',timings),('baseline_continuous',baselines),('baseline_categorical',cells)]:
 pd.DataFrame(rows).to_csv(P/f'audit/{name}.csv',index=False,encoding='utf-8-sig')
(P/'audit/data_sources.json').write_text(json.dumps(sources,ensure_ascii=False,indent=2),encoding='utf-8')
(P/'audit/extraction_manifest.json').write_text(json.dumps({'timestamp':datetime.now().astimezone().isoformat(),'checks':len(checks),'all_passed':all(c['passed'] for c in checks),'python':sys.version,'pandas':pd.__version__,'pyreadstat':pyreadstat.__version__},indent=2),encoding='utf-8')
print('All extraction checks passed:',len(checks))
