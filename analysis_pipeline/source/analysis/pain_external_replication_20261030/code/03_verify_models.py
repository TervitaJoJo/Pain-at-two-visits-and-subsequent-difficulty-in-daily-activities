from pathlib import Path
import json,sys
import pandas as pd,numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
sys.stdout.reconfigure(encoding='utf-8')
P=Path(__file__).resolve().parents[1]
rows=[]
for c in ['CHARLS','SHARE']:
 p=P/f'data_private/{c}_first_imputation_intervals.csv'
 if not p.exists():continue
 d=pd.read_csv(p,dtype={'exposure':str})
 form='y ~ C(exposure) + age + C(female) + C(education) + C(smoking) + C(arthritis) + C(cvd) + cesd + C(interval)'
 if c=='SHARE':form+=' + C(country)'
 m=smf.glm(form,data=d,family=sm.families.Binomial(link=sm.families.links.CLogLog()),offset=np.log(d.years)).fit(cov_type='cluster',cov_kwds={'groups':d.pid,'use_correction':True},maxiter=100,tol=1e-10)
 X=m.model.exog;mu=np.asarray(m.fittedvalues);weights=m.model.family.weights(mu)
 Ainv=np.linalg.inv(X.T@(weights[:,None]*X))
 scores=X*(((np.asarray(d.y)-mu)/m.model.family.link.deriv(mu)/m.model.family.variance(mu))[:,None])
 grouped=pd.DataFrame(scores).groupby(d.pid.to_numpy()).sum().to_numpy();G=len(grouped);N,K=X.shape
 V=Ainv@(grouped.T@grouped)@Ainv*(G/(G-1))*((N-1)/(N-K))
 se_fisher=pd.Series(np.sqrt(np.diag(V)),index=m.params.index)
 rr=pd.read_csv(P/f'results/{c}_first_imputation_coefficients.csv').set_index('term')
 for g in ['01','10','11']:
  key=f'C(exposure)[T.{g}]';ref=rr.loc[f'exposure{g}'];err=abs(m.params[key]-ref.estimate);ratio=se_fisher[key]/ref.se
  ok=err<1e-4 and abs(ratio-1)<.0001
  rows.append(dict(cohort=c,contrast=f'{g}_vs_00',R_coefficient=float(ref.estimate),Python_coefficient=float(m.params[key]),absolute_coefficient_difference=float(err),Python_R_SE_ratio=float(ratio),passed=bool(ok)))
  assert ok,(c,g,err,ratio)
 print(c,'independent Python GLM coefficients and cluster SE check passed',flush=True)
p=P/'audit/independent_model_verification.json';p.write_text(json.dumps({'checks':rows,'note':'Same completed dataset, independent Python GLM solver and manual participant-cluster sandwich using expected Fisher-information bread, matching R GLM HC1. Coefficients and SE ratios must agree to 0.0001.'},indent=2),encoding='utf-8')

