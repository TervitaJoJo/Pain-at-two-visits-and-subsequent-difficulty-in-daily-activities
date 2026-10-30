from pathlib import Path
import json, hashlib, platform
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import minimize_scalar
from scipy.stats import t,norm,chi2
P=Path(__file__).resolve().parents[1]
d=pd.read_csv(P/'results/meta_inputs.csv'); out=pd.read_csv(P/'results/meta_estimates.csv')
checks=[]
for row in out.to_dict('records'):
 a=d[d.cohort!=row['analysis'][5:]] if row['analysis'].startswith('omit_') else d
 y=a.logHR.to_numpy(); v=a.se.to_numpy()**2;k=len(y)
 def objective(tau):
  w=1/(v+tau); mu=np.sum(w*y)/sum(w)
  return .5*(np.log(v+tau).sum()+np.log(w.sum())+np.sum(w*(y-mu)**2))
 opt=minimize_scalar(objective,bounds=(0,2),method='bounded',options={'xatol':1e-13})
 tau=0. if objective(0)<=opt.fun else opt.x
 if row['method']=='FE':tau=0.
 w=1/(v+tau);mu=np.sum(w*y)/sum(w);q=np.sum(w*(y-mu)**2)/(k-1)
 scale=max(1,q) if row['test']=='adhoc' else q if row['test']=='knha' else 1.
 se=np.sqrt(scale/sum(w));pval=2*t.sf(abs(mu/se),k-1) if row['test'] in ['adhoc','knha'] else 2*norm.sf(abs(mu/se))
 crit=t.ppf(.975,k-1) if row['test'] in ['adhoc','knha'] else norm.ppf(.975)
 got=dict(logHR=mu,se=se,tau2=tau,HR=np.exp(mu),lower=np.exp(mu-crit*se),upper=np.exp(mu+crit*se),p=pval)
 deltas={key:float(abs(value-row[key])) for key,value in got.items()}
 passed=all(np.isclose(value,row[key],rtol=1e-4,atol=1e-8) for key,value in got.items())
 assert passed,(row['analysis'],got,deltas)
 checks.append({'analysis':row['analysis'],'passed':bool(passed),'absolute_differences':deltas})
w0=1/d.se.to_numpy()**2;y=d.logHR.to_numpy();mu0=np.sum(w0*y)/sum(w0);Q=np.sum(w0*(y-mu0)**2)
lock=json.loads((P/'audit/specification_lock.json').read_text()); current=hashlib.sha256((P/'audit/analysis_specification.md').read_text(encoding='utf-8').encode('utf-8')).hexdigest();assert lock['sha256']==current
assert sum(d.persons)==6998 and sum(d.events)==1232 and abs(d.weight_percent.sum()-100)<1e-10
res={'method':'Independent scipy bounded REML optimization and explicit modified-HK/HK/Wald calculations','python':platform.python_version(),'scipy':scipy.__version__,'checks':checks,'Q_based_I2_percent':float(max(0,(Q-3)/Q)*100),'REML_based_I2_in_main_output':float(out.iloc[0].I2),'plan_hash_unchanged':True,'plan_hash_normalization':'UTF-8 text with LF newlines, matching the initial specification lock','source_hashes':{str(Path(f)):hashlib.sha256((P.parents[2]/f).read_bytes()).hexdigest() for f in d.source}}
(P/'audit/independent_meta_verification.json').write_text(json.dumps(res,indent=2),encoding='utf-8')
print('Verified',len(checks),'model summaries; exact MI SEs and cohort counts retained.')
print('Q-based I2:',res['Q_based_I2_percent'],'; REML-based I2:',res['REML_based_I2_in_main_output'])
