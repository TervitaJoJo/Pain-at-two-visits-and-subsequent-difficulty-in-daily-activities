# Adapted from the existing Vivid empirical.forest working copy.
# CI segments, cap style, cohort markers, log scale and numerical columns are preserved.
# Two vertical panels separate distinct estimands; no pooled effect is implied.
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FixedFormatter,NullLocator
from PIL import Image
from _utils import plot_utils as pu
from figure_data import D
pu.setup_style();PALETTE,COLORS=pu.PALETTE,pu.COLORS
plt.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
fig=plt.figure(figsize=(170/25.4,120/25.4))
specs=[('A  First-visit severity adjustment',[(c,'severity_adjusted_MI','10_vs_11_severity',c+' | No longer reported vs both') for c in ['HRS','ELSA']],[.425,.65,.29,.20]),('B  Follow-up intervals: both visits vs neither',[(c,'interval'+str(i)+'_diagnostic_MI','11_vs_00',c+' | Interval '+str(i)) for c in ['HRS','ELSA'] for i in [1,2]],[.425,.12,.29,.36])]
for title,rows,rect in specs:
 ax=fig.add_axes(rect);tr=ax.get_yaxis_transform()
 for j,(c,variant,ct,label) in enumerate(rows):
  r=D[(D.cohort==c)&(D.variant==variant)&(D.contrast==ct)].iloc[0]
  if j%2==0:ax.axhspan(j-.44,j+.44,color=COLORS['bg_box'],zorder=0)
  ax.plot([r.lower,r.upper],[j,j],color=COLORS['text'],lw=1.2,solid_capstyle='round')
  for end in [r.lower,r.upper]:ax.plot([end,end],[j-.09,j+.09],color=COLORS['text'],lw=1)
  ax.plot(r.estimate,j,'o',markersize=5.5,markerfacecolor=PALETTE[0] if c=='HRS' else 'white',markeredgecolor=PALETTE[0] if c=='HRS' else PALETTE[1],markeredgewidth=1.3)
  ax.text(-.06,j,label,transform=tr,ha='right',va='center',fontsize=9,clip_on=False)
  ax.text(1.08,j,f'{r.estimate:.2f} ({r.lower:.2f}–{r.upper:.2f})',transform=tr,ha='left',va='center',fontsize=9,clip_on=False)
 ax.axvline(1,ls='--',lw=1.1,color=COLORS['ref_line'],alpha=.7)
 ax.set_xscale('log');ax.set_xlim(.3,5.2);ax.set_ylim(len(rows)-.4,-.6);ax.set_yticks([])
 ax.xaxis.set_major_locator(FixedLocator([.5,1,2,4]));ax.xaxis.set_major_formatter(FixedFormatter(['0.5','1','2','4']));ax.xaxis.set_minor_locator(NullLocator())
 ax.tick_params(axis='x',labelsize=8.5)
 for spine in ['top','right','left']:ax.spines[spine].set_visible(False)
 ax.grid(axis='x',alpha=.15,ls='--')
 ax.set_xlabel('Adjusted hazard ratio',fontsize=9)
 fig.text(.02,rect[1]+rect[3]+.05,title,fontsize=9.5,fontweight='bold',va='center')
fig.text(.738,.955,'Estimate (95% CI)',fontsize=9.5,fontweight='bold',va='center')
out=Path(__file__).resolve().parents[1]/'supplementary_figures';base='Figure_S1_secondary_diagnostics'
fig.savefig(out/(base+'.pdf'),facecolor='white');fig.savefig(out/(base+'.png'),dpi=600,facecolor='white')
with Image.open(out/(base+'.png')) as im:im.convert('RGB').save(out/(base+'.tif'),compression='tiff_lzw',dpi=(600,600))
plt.close(fig)
print('Two-panel diagnostics exported')