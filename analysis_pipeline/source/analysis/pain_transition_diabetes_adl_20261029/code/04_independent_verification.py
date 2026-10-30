from pathlib import Path
from datetime import datetime
import pandas as pd, numpy as np, json, hashlib, sys
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
P=Path(__file__).resolve().parents[1]
checks=[]
def check(cohort,name,ok,detail=''):
 checks.append(dict(cohort=cohort,check=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(f'{cohort}: {name}: {detail}')
for cohort in ['HRS','ELSA']:
 d=pd.read_csv(P/f'data_private/{cohort}_longitudinal.csv',dtype={'exposure':str})
 rows=[]
 for _,person in d.iterrows():
  previous=0.
  for j in [1,2]:
   y=person[f'y{j}']
   if pd.isna(y):break
   row=person.to_dict();row.update(y=y,interval=j,years=person[f't{j}']-previous)
   rows.append(row)
   if y==1:break
   previous=person[f't{j}']
 long=pd.DataFrame(rows)
 expected_n,expected_events=(2650,510) if cohort=='HRS' else (544,99)
 check(cohort,'independent_interval_person_count',long.pid.nunique()==expected_n)
 check(cohort,'independent_interval_event_count',long.y.sum()==expected_events)
 check(cohort,'at_most_one_first_event',long.groupby('pid').y.sum().max()<=1)
 check(cohort,'positive_time',long.years.gt(0).all())
 cc=long.dropna(subset=['age','female','education','smoking','arthritis','cvd','cesd']).copy()
 formula='y ~ C(exposure, Treatment(reference="00")) + age + C(female) + C(education) + C(smoking) + C(arthritis) + C(cvd) + cesd + C(entry) + C(interval)'
 fit=smf.glm(formula,data=cc,family=sm.families.Binomial(link=sm.families.links.CLogLog()),offset=np.log(cc.years)).fit(maxiter=100,tol=1e-10,cov_type='cluster',cov_kwds={'groups':cc.pid,'use_correction':True})
 # Match R sandwich::bread.glm's expected-information bread. statsmodels' default
 # robust GLM uses an observed-Hessian bread for this noncanonical link.
 X=fit.model.exog; mu=fit.fittedvalues.to_numpy(); dm=fit.model.family.link.inverse_deriv(fit.predict(which="linear"))
 var=mu*(1-mu); bread=np.linalg.inv(X.T@((dm*dm/var)[:,None]*X))
 scores=X*((cc.y.to_numpy()-mu)*dm/var)[:,None]
 grouped=pd.DataFrame(scores).groupby(cc.pid.to_numpy()).sum().to_numpy()
 N,K=X.shape;G=cc.pid.nunique()
 matched_cov=bread@(grouped.T@grouped)@bread*(G/(G-1))*((N-1)/(N-K))
 contrasts=pd.read_csv(P/f'results/{cohort}_contrasts.csv')
 # The two GLM solvers use different stopping rules; require <0.01% relative effect difference.
 for pair in [('11','00'),('01','00'),('10','00'),('10','11')]:
  name='_vs_'.join(pair);v=np.zeros(len(fit.params))
  for sign,g in zip([1,-1],pair):
   if g!='00':v[list(fit.params.index).index(f'C(exposure, Treatment(reference="00"))[T.{g}]')]=sign
  est=v @ fit.params;se=np.sqrt(v@matched_cov@v)
  ref=contrasts[(contrasts.variant=='complete_case')&(contrasts.contrast==name)].iloc[0]
  check(cohort,f'python_R_CC_coefficient_{name}',abs(est-ref.log_estimate)<1e-4,abs(est-ref.log_estimate))
  check(cohort,f'python_R_CC_cluster_SE_{name}',abs(se-ref.se)/ref.se<1e-4,abs(se-ref.se)/ref.se)
 unpooled=pd.read_csv(P/f'audit/{cohort}_primary_unpooled_contrasts.csv')
 for name,z in unpooled.groupby('contrast'):
  ref=contrasts[(contrasts.variant=='primary_MI')&(contrasts.contrast==name)].iloc[0]
  m=len(z);q=z.Q.mean();b=z.Q.var(ddof=1);u=z.U.mean();total=u+(1+1/m)*b
  lam=max((1+1/m)*b/total,1e-4);dfcom=z.n.min()-z.k.max();dfold=(m-1)/lam**2;dfobs=(dfcom+1)/(dfcom+3)*dfcom*(1-lam);df=1/(1/dfold+1/dfobs)
  ci=np.exp(q+np.array([-1,1])*stats.t.ppf(.975,df)*np.sqrt(total))
  check(cohort,f'independent_Rubin_pool_{name}',abs(q-ref.log_estimate)<1e-10 and abs(np.sqrt(total)-ref.se)<1e-10)
  check(cohort,f'independent_Rubin_CI_{name}',np.max(abs(ci-np.array([ref.lower,ref.upper])))<1e-8)
 diag=pd.read_csv(P/f'audit/{cohort}_model_diagnostics.csv')
 check(cohort,'all_glm_converged',diag.converged.all())
 check(cohort,'robust_covariance_positive',diag.min_cov_eigen.min()>0)
 check(cohort,'all_primary_models_same_n_events',diag.loc[diag.variant.eq('primary_MI'),['persons','events']].drop_duplicates().shape[0]==1)
 pri=contrasts[contrasts.variant.eq('primary_MI')]
 check(cohort,'all_four_prespecified_contrasts_present',len(pri)==4)
 check(cohort,'MI_Monte_Carlo_error_small',np.max(pri.mcse_log/pri.se)<.05,np.max(pri.mcse_log/pri.se))
 for sample in ['longitudinal','cross_sectional']:
  conv=pd.read_csv(P/f'audit/{cohort}_{sample}_mice_convergence.csv');end=conv[conv['.it'].eq(20)].psrf.dropna()
  check(cohort,f'{sample}_imputation_psrf',len(end)>0 and end.max()<1.1,end.max())
 # Only aggregates are written into verification outputs.
 print(cohort,'independent Python/R agreement passed; pooled contrasts verified',flush=True)
plan=P/'audit/frozen_analysis_plan_v1.md';lock=json.loads((P/'audit/plan_lock.json').read_text(encoding='utf-8'))
check('both','frozen_plan_hash_unchanged',hashlib.sha256(plan.read_bytes()).hexdigest()==lock['sha256'])
pd.DataFrame(checks).to_csv(P/'audit/independent_verification.csv',index=False,encoding='utf-8-sig')
(P/'audit/verification_summary.json').write_text(json.dumps({'timestamp':datetime.now().astimezone().isoformat(),'checks':len(checks),'passed':all(x['passed'] for x in checks),'Python':sys.version,'statsmodels':sm.__version__},indent=2),encoding='utf-8')
print('Independent verification:',len(checks),'checks passed')
