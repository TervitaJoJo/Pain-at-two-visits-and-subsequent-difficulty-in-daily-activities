from pathlib import Path
import pandas as pd
import matplotlib as mpl
# Adapted from the complete vivid advanced.dumbbell recipe.
# Native width 180 mm; displayed at 180 mm, so 9 pt remains 9 pt.
P=Path(__file__).resolve().parents[1]
PALETTE=['#606C38','#DDA15E','#BC6C25','#A68A64','#936639','#7F4F24','#4A5A2B']
COLORS={'grid':'#cccccc'}
mpl.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42,'ps.fonttype':42})
import matplotlib.pyplot as plt
import numpy as np

metrics=['HRS','ELSA','SHARE']
d=pd.read_csv(P/'results/standardized_net_risks.csv');d=d[d.variant=='primary_MI'].set_index(['cohort','estimand'])
before=[100*d.loc[(c,'risk_00'),'estimate'] for c in metrics]
after=[100*d.loc[(c,'risk_11'),'estimate'] for c in metrics]

# ── 自适应高度
_fig_h = 3.15
fig, ax = plt.subplots(figsize=(180/25.4, _fig_h))
y = np.arange(len(metrics))

# Subtle grid
ax.grid(axis='x', alpha=0.15, linestyle='-', color=COLORS['grid'])
ax.set_axisbelow(True)

# Simple connector lines + arrow heads
for i in range(len(metrics)):
    delta = after[i] - before[i]
    pct_change = delta / before[i] * 100
    # Simple connector line
    ax.plot([before[i], after[i]], [y[i], y[i]],
            color=PALETTE[2], linewidth=1.6, alpha=.7, solid_capstyle='round')
    # Arrow head at the "after" end
    ax.annotate('', xy=(after[i], y[i]), xytext=(after[i] - 1.8, y[i]),
                arrowprops=dict(arrowstyle='->', color=PALETTE[0], lw=2.0))

    # Before / After dots
    ax.scatter(before[i], y[i], color=PALETTE[3], s=80, zorder=3,
               edgecolors='white', linewidths=1.0, label='No pain at either visit (00)' if i == 0 else '')
    ax.scatter(after[i], y[i], color=PALETTE[0], s=80, zorder=3,
               edgecolors='white', linewidths=1.0, label='Pain at both visits (11)' if i == 0 else '')

    # Model-based 95% CIs and associated risk difference, not a treatment effect.
    for tag,color in [('risk_00',PALETTE[3]),('risk_11',PALETTE[0])]:
        r=d.loc[(metrics[i],tag)]
        ax.errorbar(r['estimate']*100,y[i],xerr=[[100*(r['estimate']-r['lower'])],[100*(r['upper']-r['estimate'])]],fmt='none',ecolor=color,elinewidth=1.15,capsize=3,zorder=2)
    r=d.loc[(metrics[i],'RD_11_vs_00')]
    ax.text(37,y[i],f"{r['estimate']*100:+.1f} ({r['lower']*100:+.1f}, {r['upper']*100:+.1f})",va='center',fontsize=9,color='#333333')

ax.set_yticks(y)
ax.set_yticklabels(metrics, fontsize=10)
ax.set_xlabel('Model-standardized 4-year net ADL risk (%)', fontsize=10)
ax.set_xlim(7,36)
ax.set_xticks([10,15,20,25,30,35])
ax.set_ylim(-.6,2.55)
ax.text(37,-.62,'Risk difference (95% CI), pp',fontsize=9,va='bottom')
ax.legend(loc='upper left', bbox_to_anchor=(0,-.22), frameon=False, labelspacing=0.35, handlelength=1.6, fontsize=9,
          fancybox=True, shadow=False)
ax.invert_yaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.subplots_adjust(left=.115,right=.67,top=.84,bottom=.30)
fig.text(.115,.94,'Two-visit pain patterns and subsequent daily-activity difficulty',fontsize=11,weight='bold')
for ext in ['png','tif','pdf']:
    fig.savefig(P/'figures'/f'standardized_absolute_risks.{ext}',dpi=600,facecolor='white')
plt.close(fig)
