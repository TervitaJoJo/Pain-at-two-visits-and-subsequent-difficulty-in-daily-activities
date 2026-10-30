import numpy as np, matplotlib.pyplot as plt
from _utils.plot_utils import setup_style, save_fig, PALETTE, COLORS, _lighten
setup_style()

from pathlib import Path
import pandas as pd
from _utils import plot_utils as pu
P = Path(__file__).resolve().parents[1]
PALETTE, COLORS = pu.PALETTE, pu.COLORS
plt.rcParams.update({'font.family': 'Arial', 'font.size': 9, 'pdf.fonttype': 42})
d = pd.read_csv(P / 'results/meta_inputs.csv')
r = pd.read_csv(P / 'results/meta_estimates.csv').set_index('analysis').loc['all_four_modified_HK']
labels = d.cohort.replace({'CHARLS': 'CHARLS*'}).tolist()
coefs, ci_lo, ci_hi = d.HR.to_numpy(), d.lower.to_numpy(), d.upper.to_numpy()
weights = d.weight_percent.to_numpy()
significant = (ci_lo > 1) | (ci_hi < 1)
pooled_coef, pooled_lo, pooled_hi = r.HR, r.lower, r.upper
# Native 7.2-inch width = 182.9 mm; at 180 mm, 9 pt text appears at 8.86 pt.
# Rows/CI caps, weight-scaled open/filled circles, aligned columns and diamond
# retained from empirical.forest. Null line becomes HR=1 on a log axis.
# ★ 自适应高度：变量多时不会过度拉伸
_fig_h = 4.4
fig, ax = plt.subplots(figsize=(7.2, _fig_h))
y = np.arange(len(labels)) * 1.4

# 交替行阴影
for i in range(len(labels)):
    if i % 2 == 0:
        ax.axhspan(y[i] - 0.55, y[i] + 0.55, color=COLORS['bg_box'], zorder=0)

# 零线参考
ax.axvline(x=1, color=COLORS['down'], linestyle='--', linewidth=1.5, alpha=0.7, zorder=1)

for i, (c, lo, hi, lbl, sig, wt) in enumerate(zip(coefs, ci_lo, ci_hi, labels, significant, weights)):
    # CI 线 + 端点帽
    ax.plot([lo, hi], [y[i], y[i]], color=COLORS['text'], linewidth=1.2, zorder=2, solid_capstyle='round')
    ax.plot([lo, lo], [y[i]-0.12, y[i]+0.12], color=COLORS['text'], linewidth=1.0)
    ax.plot([hi, hi], [y[i]-0.12, y[i]+0.12], color=COLORS['text'], linewidth=1.0)

    # 权重比例圆点（显著=实心，不显著=空心）
    fc = PALETTE[0] if sig else 'white'
    ec = PALETTE[0] if sig else COLORS['ref_line']
    ms = 10 * np.sqrt(wt / max(weights))
    ax.plot(c, y[i], 'o', color=fc, markersize=ms, markeredgecolor=ec,
            markeredgewidth=1.3, zorder=3)

    # 左侧标签
    ax.text(0.73, y[i], lbl, ha='right', va='center', fontsize=9, color=COLORS['text'])

    # 右侧数值列：系数 [95% CI] + 显著性★
    sig_mark = ''
    ax.text(3.45, y[i], f'{c:.2f} ({lo:.2f}-{hi:.2f})',
            ha='left', va='center', fontsize=8, fontfamily='Arial',
            color=COLORS['text'],
            fontweight='normal')

    # 权重列
    ax.text(8.1, y[i], f'{wt:.1f}%', ha='right', va='center',
            fontsize=8, color=COLORS['text'])

# Aligned column headings.
for x, label, align in [(0.73, 'Cohort', 'right'), (3.45, 'Adjusted HR (95% CI)', 'left'), (8.1, 'Weight', 'right')]:
    ax.text(x, -1.05, label, ha=align, va='center', fontsize=9, fontweight='bold', color=COLORS['text'])
# 菱形 pooled estimate
diamond_y = max(y) + 1.6
dh = 0.35
diamond_x = [pooled_lo, pooled_coef, pooled_hi, pooled_coef]
diamond_yy = [diamond_y, diamond_y + dh, diamond_y, diamond_y - dh]
ax.fill(diamond_x, diamond_yy, color=PALETTE[1], alpha=0.85, zorder=4)
ax.plot(diamond_x + [diamond_x[0]], diamond_yy + [diamond_yy[0]],
        color=COLORS['text'], linewidth=0.8, zorder=5)
ax.text(0.73, diamond_y, 'Pooled', ha='right', va='center',
        fontsize=9.5, fontweight='bold', color=PALETTE[1])
ax.text(3.45, diamond_y, f'{pooled_coef:.2f} ({pooled_lo:.2f}-{pooled_hi:.2f})',
        ha='left', va='center', fontsize=8, fontfamily='Arial',
        fontweight='bold', color=PALETTE[1])

# Model and heterogeneity summary from the fitted model.
i2_text = f"REML + modified Hartung-Knapp: P = {r.p:.4f};  I$^2$ = {r.I2:.1f}%"
ax.text(0.70, diamond_y + 1.05, i2_text, ha='left', va='center', fontsize=8.5,
        bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor=COLORS['ref_line']),
        color=COLORS['text'])
ax.set_yticks([])
ax.set_xlabel('Adjusted hazard ratio (log scale)', fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_ylim(-1.7, diamond_y + 1.7)
ax.set_xscale('log')
ax.set_xlim(0.68, 8.45)
ax.set_xticks([0.8, 1, 1.5, 2, 3])
ax.set_xticklabels(['0.8', '1', '1.5', '2', '3'])
ax.minorticks_off()
ax.invert_yaxis()
ax.grid(axis='x', alpha=0.15, linestyle='--')
fig.suptitle('Pain at both visits and incident ADL difficulty', x=0.53, y=0.98, fontsize=11, fontweight='bold')
fig.text(0.03, 0.025, '*CHARLS: adapted pain recall. Open circles: cohort CI includes 1.', fontsize=8)
fig.tight_layout(rect=[0, 0.06, 1, 0.95])
for ext in ['png', 'pdf', 'tif']:
    kw = {'pil_kwargs': {'compression': 'tiff_lzw'}} if ext == 'tif' else {}
    fig.savefig(P / 'figures' / f'four_cohort_meta_forest.{ext}', dpi=600, facecolor='white', **kw)
plt.close(fig)

