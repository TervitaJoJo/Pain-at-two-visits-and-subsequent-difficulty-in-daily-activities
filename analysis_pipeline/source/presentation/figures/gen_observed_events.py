from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
from _utils import plot_utils as pu
pu.setup_style()
PALETTE,COLORS=pu.PALETTE,pu.COLORS
_lighten=pu._lighten
P=Path(__file__).resolve().parents[1]
A=Path(__file__).resolve().parents[2]/'analysis/pain_transition_diabetes_adl_20261029'
source=pd.read_csv(A/'results/table2_followup_by_group.csv',dtype={'group':str})
codes=['00','01','10','11']
categories=['Neither visit\n(00)','Newly reported\n(01)','No longer reported\n(10)','Both visits\n(11)']
groups={c:source[source.cohort.eq(c)].set_index('group').loc[codes,'event_fraction_among_observed_percent'].to_numpy() for c in ['HRS','ELSA']}
# Native 170 mm and placement 170 mm: 9 pt remains 9 pt on the manuscript page.
# Retain recipe's light fills, colored outlines and offset shadows.
# No best-method star or pooled mean: cohorts are not methods to be ranked.
# No error bars: observed fractions have unequal observation time and are shown descriptively.
plt.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
fig,ax=plt.subplots(figsize=(170/25.4,112/25.4))
fig.subplots_adjust(left=.12,right=.985,bottom=.22,top=.88)
x=np.arange(len(categories));n=len(groups);width=.33
ax.axhspan(0,39,alpha=.03,color=PALETTE[0],zorder=0)
export=[]
for i,(name,vals) in enumerate(groups.items()):
 offset=(i-n/2+.5)*width
 ax.bar(x+offset+.02,vals,width,color='#cccccc',alpha=.08,zorder=1)
 bars=ax.bar(x+offset,vals,width,color=_lighten(PALETTE[i],.4),edgecolor=PALETTE[i],linewidth=1.2,label=name,zorder=2)
 rows=source[source.cohort.eq(name)].set_index('group').loc[codes]
 for j,(bar,v) in enumerate(zip(bars,vals)):
  r=rows.iloc[j]
  assert np.isclose(v,100*r.events/r.with_followup)
  ax.text(bar.get_x()+bar.get_width()/2,bar.get_height()+.75,f'{v:.1f}%',ha='center',va='bottom',fontsize=9,color=COLORS['text'])
  ax.text(bar.get_x()+bar.get_width()/2,bar.get_height()/2,f'{int(r.events)} /\n{int(r.with_followup):,}',ha='center',va='center',fontsize=8,color=COLORS['text'],linespacing=1.2)
  export.append(dict(cohort=name,pattern=codes[j],events=int(r.events),participants_with_followup=int(r.with_followup),observed_percent=float(v)))
ax.set_xticks(x);ax.set_xticklabels(categories,fontsize=9,linespacing=1.4)
ax.set_ylabel('Participants with an observed ADL event (%)',fontsize=9)
ax.legend(frameon=False,handlelength=1.6,fontsize=9,loc='upper left',ncols=2,bbox_to_anchor=(0,1.17))
ax.set_ylim(0,39);ax.set_yticks([0,10,20,30]);ax.grid(axis='y',alpha=.12,linestyle='--',color=COLORS['grid']);ax.set_axisbelow(True)
ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
fig.text(.12,.035,'Within bars: first events / participants with follow-up',ha='left',fontsize=8.5,color=COLORS['text'])
out=P/'main_figures';base=out/'Figure_2_observed_events'
fig.savefig(base.with_suffix('.pdf'),facecolor='white')
fig.savefig(base.with_suffix('.png'),dpi=600,facecolor='white')
with Image.open(base.with_suffix('.png')) as im:im.convert('RGB').save(base.with_suffix('.tif'),compression='tiff_lzw',dpi=(600,600))
plt.close(fig)
(P/'qa').mkdir(exist_ok=True)
pd.DataFrame(export).to_csv(P/'qa'/'figure2_plotted_data.csv',index=False)
print('Rendered Figure 2; all eight bars verified against source event counts.')