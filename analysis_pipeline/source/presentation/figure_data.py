from pathlib import Path
import pandas as pd
A=Path(__file__).resolve().parents[1]/'analysis/pain_transition_diabetes_adl_20261029'
D=pd.concat([pd.read_csv(A/'results/all_contrasts.csv')]+[pd.read_csv(A/f'revision_20261029/results/{c}_revision_contrasts.csv') for c in ['HRS','ELSA']],ignore_index=True)
def forest_data(kind):
 rows=[];headers=[];pos=0
 for c in ['HRS','ELSA']:
  headers.append({'y':pos-1,'label':c+('  |  United States' if c=='HRS' else '  |  England')})
  if kind=='primary':
   specs=[('primary_MI',k,l) for k,l in [('11_vs_00','Both visits vs neither'),('01_vs_00','Newly reported vs neither'),('10_vs_00','No longer reported vs neither'),('10_vs_11','No longer reported vs both')]]
  elif kind=='sensitivity':
   specs=[(k,'11_vs_00',l) for k,l in [('complete_case','Complete cases'),('function_adjusted_MI','Function at first visit'),('landmark_function_adjusted_MI','Function at landmark'),('treatment_adjusted_MI','Diabetes treatment')]]
   if c=='HRS':specs.append(('diabetes_context_MI','11_vs_00','Treatment + BMI + duration'))
   specs += [(k,'11_vs_00',l) for k,l in [('proxy_censored_original_MI','Censored before proxy'),('observation_weighted_MI','Observation weighting'),('nonlinear_covariates_MI','Nonlinear age and CES-D8')]]
  elif kind=='diagnostics':
   specs=[('severity_adjusted_MI','10_vs_11_severity','Severity-adjusted: 10 vs 11'),('interval1_diagnostic_MI','11_vs_00','First interval: 11 vs 00'),('interval2_diagnostic_MI','11_vs_00','Second interval: 11 vs 00')]
  else:specs=[('cross_sectional_MI','pain_vs_no_pain','Multiple imputation'),('cross_sectional_CC','pain_vs_no_pain','Complete cases')]
  for k,contrast,label in specs:
   d=D[(D.cohort==c)&(D.variant==k)&(D.contrast==contrast)]
   assert len(d)==1,(c,k,contrast)
   r=d.iloc[0].to_dict();r.update(label=label,y=pos);rows.append(r);pos+=1.4
  pos+=1.65
 if kind=='primary':name='Figure_3_primary';folder='main_figures';xl='Adjusted hazard ratio';lims=(.4,4.2);ticks=[.5,1,2,4]
 elif kind=='sensitivity':name='Figure_S3_sensitivity';folder='supplementary_figures';xl='Adjusted hazard ratio';lims=(.5,4.2);ticks=[.5,1,2,4]
 elif kind=='diagnostics':name='Figure_S1_secondary_diagnostics';folder='supplementary_figures';xl='Adjusted hazard ratio';lims=(.3,5.2);ticks=[.5,1,2,4]
 else:name='Figure_S2_cross_sectional';folder='supplementary_figures';xl='Adjusted prevalence ratio';lims=(.8,3.4);ticks=[1,1.5,2,3]
 return rows,headers,name,folder,xl,lims,ticks
