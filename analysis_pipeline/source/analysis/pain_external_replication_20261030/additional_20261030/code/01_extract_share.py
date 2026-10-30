import os
DATA_ROOT=__import__("pathlib").Path(os.environ["COHORT_DATA_ROOT"]).resolve()
from pathlib import Path
import json
import numpy as np
import pandas as pd
import pyreadstat
P=Path(__file__).resolve().parents[1]; ROOT=P.parents[2]
def mod(w, s, cols):
    path=next((DATA_ROOT/'SHARE').rglob(f'sharew{w}_rel9-0-0_{s}.dta'))
    d,_=pyreadstat.read_dta(path,usecols=['mergeid']+cols)
    assert not d.mergeid.duplicated().any()
    return d.set_index('mergeid')
def adl(d):
    x=d[[f'ph049d{i}' for i in [1,3,4,5,6]]]
    return pd.Series(np.where(x.eq(1).any(axis=1),1,np.where(x.eq(0).all(axis=1),0,np.nan)),index=d.index)
phcols=['country','ph006d5','ph084_','ph054_']+[f'ph049d{i}' for i in [1,3,4,5,6]]
a=mod(5,'ph',phcols+['ph009_5','ph011d6']); b=mod(6,'ph',phcols+['ph049d2','ph048d3','ph048d5'])
ac=mod(5,'cv_r',['age_int','int_year','int_month']);bc=mod(6,'cv_r',['int_year','int_month'])
idx=a.index.sort_values();a=a.reindex(idx);b=b.reindex(idx);ac=ac.reindex(idx);bc=bc.reindex(idx)
am=ac.int_year*12+ac.int_month;bm=bc.int_year*12+bc.int_month
eligible=ac.age_int.ge(50)&ac.int_month.between(1,12)&bc.int_month.between(1,12)&a.ph006d5.eq(1)&a.ph054_.isin([1,2])&b.ph054_.isin([1,2])&(bm-am).between(12,36)&adl(a).eq(0)&adl(b).eq(0)&a.ph084_.isin([1,5])&b.ph084_.isin([1,5])
ids=idx[eligible];d=pd.read_csv(P.parent/'data_private/SHARE_primary_longitudinal.csv',dtype={'exposure':str})
assert len(ids)==len(d)==4215
for name,vals in [('age',ac.loc[ids,'age_int']),('country',a.loc[ids,'country']),('landmark_month',bm.loc[ids])]:
    np.testing.assert_allclose(d[name].to_numpy(dtype=float),vals.to_numpy(dtype=float))
exp=a.loc[ids,'ph084_'].map({1:'1',5:'0'})+b.loc[ids,'ph084_'].map({1:'1',5:'0'})
assert exp.tolist()==d.exposure.tolist()
d['diabetes_medication']=a.loc[ids,'ph011d6'].where(a.loc[ids,'ph011d6'].isin([0,1])).to_numpy()
h=mod(5,'gv_health',['bmi']).reindex(ids)
d['bmi']=h.bmi.where(h.bmi.between(10,70)).to_numpy()
age_dx=a.loc[ids,'ph009_5'].where(a.loc[ids,'ph009_5'].ge(0))
diagnosis_sources=[dict(wave=5,newly_recovered=int(age_dx.notna().sum()))]
for previous_wave in [4,2,1]:
    previous=mod(previous_wave,'ph',['ph009_5']).reindex(ids).ph009_5
    valid=previous.where(previous.ge(0)&previous.le(ac.loc[ids,'age_int']))
    diagnosis_sources.append(dict(wave=previous_wave,newly_recovered=int((age_dx.isna()&valid.notna()).sum())))
    age_dx=age_dx.fillna(valid)
pd.DataFrame(diagnosis_sources).to_csv(P/'audit/SHARE_diagnosis_provenance.csv',index=False)
dur=ac.loc[ids,'age_int']-age_dx
d['diabetes_duration']=dur.where(age_dx.ge(0)&dur.ge(0)&dur.le(ac.loc[ids,'age_int']+1)).to_numpy()
f=b.loc[ids,['ph049d2','ph048d3','ph048d5']]
d['landmark_function3']=f.where(f.isin([0,1])).sum(axis=1,min_count=3).to_numpy()
rows=[];death=[]
for j,w in enumerate([7,8],1):
    ph=mod(w,'ph',['ph054_']).reindex(ids);cv=mod(w,'cv_r',['deceased','int_year','int_month']).reindex(ids)
    xt=mod(w,'xt',['xt008_','xt009_']).reindex(ids)
    d[f'proxy{j}']=ph.ph054_.map({1:0.,2:1.,3:1.}).to_numpy()
    known=cv.deceased.eq(1)|xt.xt009_.between(2000,2026)
    d[f'dead{j}']=np.where(known,1,np.where(cv.deceased.eq(0),0,np.nan))
    d[f'death_month{j}']=(12*xt.xt009_+xt.xt008_).where(xt.xt009_.between(2000,2026)&xt.xt008_.between(1,12)).to_numpy()
    assert not ((d[f'dead{j}']==1)&d[f'y{j}'].notna()).any(), 'Vital status conflict'
    death.append(dict(wave=w,cv_known=int(cv.deceased.isin([0,1]).sum()),deaths=int(known.sum()),exit_year=int(xt.xt009_.between(2000,2026).sum()),exit_full_date=int(d[f'death_month{j}'].notna().sum()),outcome_known=int(d[f'y{j}'].notna().sum()),proxy_participation=int(d[f'proxy{j}'].eq(1).sum()),unknown_proxy_with_ADL=int((d[f'proxy{j}'].isna()&d[f'y{j}'].notna()).sum())))
# Death is absorbing: no reconnection after a recorded first-wave death.
d.loc[d.dead1.eq(1),'dead2']=1
for c in ['diabetes_medication','bmi','diabetes_duration','landmark_function3']:
    x=d[c];rows.append(dict(variable=c,n=len(x),missing=int(x.isna().sum()),missing_percent=100*x.isna().mean(),minimum=x.min(),maximum=x.max()))
d.to_csv(P/'data_private/SHARE_extended.csv',index=False)
pd.DataFrame(rows).to_csv(P/'audit/SHARE_extra_missingness.csv',index=False)
pd.DataFrame(death).to_csv(P/'audit/SHARE_death_proxy_audit.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(death).to_string(index=False))
print('Exact reconstruction of frozen SHARE sample passed')
