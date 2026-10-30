from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np, matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
from PIL import Image
from _utils import plot_utils as pu
from figure_data import forest_data
pu.setup_style()
PALETTE, COLORS = pu.PALETTE, pu.COLORS
plt.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
rows, headers, basename, output_folder, xlabel, xlimits, xticks = forest_data('primary')
coefs = np.array([r['estimate'] for r in rows])
ci_lo = np.array([r['lower'] for r in rows])
ci_hi = np.array([r['upper'] for r in rows])
labels = [r['label'] for r in rows]
cohorts = [r['cohort'] for r in rows]
assert np.all((ci_lo <= coefs) & (coefs <= ci_hi))

# Publication-sized canvas; data-driven rows preserve the recipe forest structure.
fig, ax = plt.subplots(figsize=(170/25.4, 118/25.4))
fig.subplots_adjust(left=0.425, right=0.715, bottom=0.14, top=0.88)
y = np.array([r['y'] for r in rows])
transform = ax.get_yaxis_transform()

# 交替行阴影
for i in range(len(labels)):
    if i % 2 == 0:
        ax.axhspan(y[i] - 0.55, y[i] + 0.55, color=COLORS['bg_box'], zorder=0)

# 零线参考
ax.axvline(x=1, color=COLORS['ref_line'], linestyle='--', linewidth=1.1, alpha=0.7, zorder=1)

for i, (c, lo, hi, lbl, cohort) in enumerate(zip(coefs, ci_lo, ci_hi, labels, cohorts)):
    # CI 线 + 端点帽
    ax.plot([lo, hi], [y[i], y[i]], color=COLORS['text'], linewidth=1.2, zorder=2, solid_capstyle='round')
    ax.plot([lo, lo], [y[i]-0.12, y[i]+0.12], color=COLORS['text'], linewidth=1.0)
    ax.plot([hi, hi], [y[i]-0.12, y[i]+0.12], color=COLORS['text'], linewidth=1.0)

    # 权重比例圆点（显著=实心，不显著=空心）
    fc = PALETTE[0] if cohort == 'HRS' else 'white'
    ec = PALETTE[0] if cohort == 'HRS' else PALETTE[1]
    ms = 5.5
    ax.plot(c, y[i], 'o', color=fc, markersize=ms, markeredgecolor=ec,
            markeredgewidth=1.3, zorder=3)

    # 左侧标签
    ax.text(-0.06, y[i], lbl, transform=transform, ha='right', va='center', fontsize=9, color=COLORS['text'], clip_on=False)

    # Aligned numerical column; filling identifies cohort, not a P-value.
    ax.text(1.08, y[i], f'{c:.2f} ({lo:.2f}–{hi:.2f})',
            transform=transform, ha='left', va='center', fontsize=9,
            color=COLORS['text'], clip_on=False)

# No pooled diamond, weights, significance stars, or heterogeneity statistic:
# the study estimates separate cohort associations, not a meta-analysis.
for h in headers:
    ax.text(-1.39, h['y'], h['label'], transform=transform, ha='left', va='center',
            fontsize=9.5, fontweight='bold', color=COLORS['text'], clip_on=False)
fig.text(0.02, 0.953, 'Comparison / model', ha='left', va='center', fontsize=9.5, fontweight='bold')
fig.text(0.738, 0.953, 'Estimate (95% CI)', ha='left', va='center', fontsize=9.5, fontweight='bold')

ax.set_yticks([])
ax.set_xlabel(xlabel, fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_ylim(-1.55, max(y) + 0.7)
ax.set_xscale('log')
ax.set_xlim(*xlimits)
ax.xaxis.set_major_locator(FixedLocator(xticks))
ax.xaxis.set_major_formatter(FixedFormatter([str(x).rstrip('0').rstrip('.') if isinstance(x,float) else str(x) for x in xticks]))
ax.xaxis.set_minor_locator(NullLocator())
ax.tick_params(axis='x',labelsize=8.5)
ax.invert_yaxis()
ax.grid(axis='x', alpha=0.15, linestyle='--')
# Preserve physical size; no tight bbox resizing.
out = Path(__file__).resolve().parents[1] / output_folder
out.mkdir(exist_ok=True)
fig.savefig(out / (basename+'.pdf'), facecolor='white')
fig.savefig(out / (basename+'.png'), dpi=600, facecolor='white')
with Image.open(out / (basename+'.png')) as im:
    im.convert('RGB').save(out / (basename+'.tif'), compression='tiff_lzw', dpi=(600,600))
plt.close(fig)
print(out / basename)
