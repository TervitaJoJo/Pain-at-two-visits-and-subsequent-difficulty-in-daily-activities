from pathlib import Path
from datetime import datetime
import json,hashlib,sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
sys.stdout.reconfigure(encoding='utf-8')
P=Path(__file__).resolve().parents[1]
R=pd.concat([pd.read_csv(P/f'results/{c}_contrasts.csv') for c in ['HRS','ELSA']],ignore_index=True)
R.to_csv(P/'results/all_contrasts.csv',index=False,encoding='utf-8-sig')
G=pd.read_csv(P/'audit/group_outcomes.csv',dtype={'group':str})
D=pd.concat([pd.read_csv(P/f'audit/{c}_model_diagnostics.csv') for c in ['HRS','ELSA']])
GNAMES={'00':'两次均无疼痛','01':'第二次新报告疼痛','10':'第二次不再报告疼痛','11':'两次均报告疼痛'}
CNAMES={'11_vs_00':'两次均有 vs 两次均无（主要）','01_vs_00':'新报告 vs 两次均无','10_vs_00':'不再报告 vs 两次均无','10_vs_11':'不再报告 vs 两次均有'}
def get(c,v='primary_MI',k='11_vs_00'):return R[(R.cohort==c)&(R.variant==v)&(R.contrast==k)].iloc[0]
def ci(r):return f'{r.estimate:.2f}（{r.lower:.2f}–{r.upper:.2f}）'
def pval(v):return '<0.001' if v<.001 else f'{v:.3f}'
def markdown(headers,rows):return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in rows)
# Publication tables contain aggregates only.
for sample in ['longitudinal','cross_sectional']:
 records=[]
 for c in ['HRS','ELSA']:
  d=pd.read_csv(P/f'data_private/{c}_{sample}.csv',dtype={'exposure':str})
  groupcol='exposure' if sample=='longitudinal' else 'pain'
  for group,z in d.groupby(groupcol):
   records.append(dict(cohort=c,group=group,variable='N',level='',value=str(len(z)),missing=0))
   for v in ['age','cesd','function3']:
    vals=z[v].dropna();records.append(dict(cohort=c,group=group,variable=v,level='',value=f'{vals.mean():.2f} ({vals.std():.2f})',missing=int(z[v].isna().sum())))
   for v in ['female','education','smoking','arthritis','cvd','pain_severity']:
    known=z[v].notna().sum()
    for level in sorted(z[v].dropna().unique()):
     n=z[v].eq(level).sum();records.append(dict(cohort=c,group=group,variable=v,level=int(level),value=f'{n} ({100*n/known:.1f}%)',missing=int(z[v].isna().sum())))
 pd.DataFrame(records).to_csv(P/f'results/table1_{sample}_baseline.csv',index=False,encoding='utf-8-sig')
longgroups=G[G['sample']=='longitudinal'].copy()
longgroups['event_fraction_among_observed_percent']=100*longgroups.events/longgroups.with_followup
longgroups['no_first_followup_percent']=100*longgroups.no_first_followup/longgroups.eligible
longgroups.to_csv(P/'results/table2_followup_by_group.csv',index=False,encoding='utf-8-sig')
primary=R[R.variant=='primary_MI'].copy();primary['estimate_95CI']=primary.apply(ci,axis=1)
primary.to_csv(P/'results/table3_primary.csv',index=False,encoding='utf-8-sig')
variants=['primary_MI','complete_case','function_adjusted_MI','observation_weighted_MI','ADL_or_death_nominal_MI','nonlinear_covariates_MI']
R[(R.contrast=='11_vs_00')&R.variant.isin(variants)].to_csv(P/'results/table4_primary_sensitivity.csv',index=False,encoding='utf-8-sig')
R[R.variant.str.startswith('cross_sectional')].to_csv(P/'results/table5_cross_sectional.csv',index=False,encoding='utf-8-sig')
