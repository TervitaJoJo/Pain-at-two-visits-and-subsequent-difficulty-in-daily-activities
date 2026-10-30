from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar
from scipy.special import expit
P=Path(__file__).resolve().parents[1];checks=[]
def check(name,ok,detail=''):
    checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError((name,detail))
def pool(z):
    m=len(z);q=z.Q.mean();u=z.U.mean();n=z.n.min();k=z.k.max()
    if m==1:return q,u,n-k
    b=z.Q.var(ddof=1);total=u+(1+1/m)*b;lam=(1+1/m)*b/total
    dc=n-k;tmp=(1-lam)*(1+dc)*dc;df=(m-1)*tmp/((dc+3)*(m-1)+lam**2*tmp)
    return q,total,df
x=pd.read_csv(P/'results/SHARE_scalar_inputs.csv');refs=pd.read_csv(P/'results/SHARE_additional_contrasts.csv')
for (v,c),z in x.groupby(['variant','contrast']):
    ref=refs[(refs.variant==v)&(refs.contrast==c)].iloc[0];q,total,df=pool(z)
    check(f'Rubin_{v}_{c}',np.allclose([q,np.sqrt(total),df],[ref.log_estimate,ref.se,ref.df],rtol=1e-9))
    ci=np.exp(q+np.array([-1,1])*stats.t.ppf(.975,df)*np.sqrt(total))
    check(f'CI_{v}_{c}',np.allclose(ci,[ref.lower,ref.upper],rtol=1e-9))
x=pd.read_csv(P/'results/standardization_scalar_inputs.csv');refs=pd.read_csv(P/'results/standardized_net_risks.csv')
for (c,v,e),z in x.groupby(['cohort','variant','estimand']):
    ref=refs[(refs.cohort==c)&(refs.variant==v)&(refs.estimand==e)].iloc[0];q,total,df=pool(z)
    ci=q+np.array([-1,1])*stats.t.ppf(.975,df)*np.sqrt(total)
    if z.scale.iloc[0]=='logit':ci=expit(ci);est=expit(q)
    else:est=q
    check(f'standardized_{c}_{v}_{e}',np.allclose([est,*ci],[ref.estimate,ref.lower,ref.upper],rtol=1e-9,atol=1e-10))
refs=pd.read_csv(P/'results/three_cohort_meta.csv');x=pd.read_csv(P/'results/three_cohort_meta_inputs.csv')
for _,r in refs.iterrows():
    z=x[(x.variant==r.variant)&(x.contrast==r.contrast)];y=z.log_estimate.to_numpy();v=z.se.to_numpy()**2;k=len(z)
    def objective(tau):
        w=1/(v+tau);mu=np.sum(w*y)/sum(w)
        return .5*(np.log(v+tau).sum()+np.log(w.sum())+np.sum(w*(y-mu)**2))
    opt=minimize_scalar(objective,bounds=(0,3),method='bounded',options={'xatol':1e-14});tau=0 if objective(0)<=opt.fun else opt.x
    w=1/(v+tau);mu=np.sum(w*y)/sum(w);scale=max(1,np.sum(w*(y-mu)**2)/(k-1));se=np.sqrt(scale/sum(w));ci=np.exp(mu+np.array([-1,1])*stats.t.ppf(.975,k-1)*se)
    check(f'REML_mHK_{r.variant}_{r.contrast}',np.allclose([mu,se,tau,*ci],[r.logHR,r.se,r.tau2,r.lower,r.upper],rtol=1e-5,atol=1e-8))
# Independent expected-Fisher clustered covariance and score check on extended context model.
f=P/'data_private/context_verification_design.csv'
if f.exists():
    d=pd.read_csv(f);r=pd.read_csv(P/'results/context_verification_coefficients.csv');X=d[r.term].to_numpy();b=r.estimate.to_numpy();eta=X@b+d.offset.to_numpy();H=np.exp(eta);mu=-np.expm1(-H);dm=H*np.exp(-H);var=mu*(1-mu)
    bread=np.linalg.inv(X.T@((dm*dm/var)[:,None]*X));score=X*((d.y.to_numpy()-mu)*dm/var)[:,None]
    group=pd.DataFrame(score).groupby(d.pid.to_numpy()).sum().to_numpy();n,k=X.shape;g=d.pid.nunique();V=bread@(group.T@group)@bread*(g/(g-1))*((n-1)/(n-k))
    check('extended_context_independent_cluster_SE',np.allclose(np.sqrt(np.diag(V)),r.se,rtol=1e-4))
    check('extended_context_score_stationarity',np.max(np.abs(score.sum(axis=0)))/n<1e-3)
for fn in ['SHARE_model_diagnostics.csv','function_model_diagnostics.csv']:
    d=pd.read_csv(P/'audit'/fn);check(fn+'_all_converged',d.converged.all());check(fn+'_no_model_warnings',not d.warnings.any())
lock=json.loads((P/'audit/specification_lock.json').read_text());h=hashlib.sha256((P/'audit/additional_analysis_specification.md').read_text(encoding='utf-8').encode()).hexdigest();check('specification_unchanged',h==lock['sha256_LF_utf8'])
(P/'audit/independent_verification.json').write_text(json.dumps(dict(checks=len(checks),all_passed=all(c['passed'] for c in checks),details=checks),indent=2),encoding='utf-8')
print(len(checks),'independent checks passed')
