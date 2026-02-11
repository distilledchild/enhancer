import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path

# This script generates a flowchart illustrating the filtering process of chromatin loops.
# It uses matplotlib to draw shapes and text, representing each step of the analysis.
# Larger default for readability while preserving layout with per-shape fitting.
BASE_FONTSIZE = 17
MIN_FONTSIZE = 11
DECISION_BASE_FONTSIZE = 21


def fit_text_fontsize(
    ax,
    x,
    y,
    text,
    fontweight,
    base_fontsize,
    max_width,
    max_height,
    linespacing=1.2,
    min_fontsize=MIN_FONTSIZE
):
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    temp = ax.text(
        x,
        y,
        text,
        ha='center',
        va='center',
        fontsize=base_fontsize,
        fontweight=fontweight,
        multialignment='center',
        linespacing=linespacing,
        alpha=0.0
    )

    fitted = float(base_fontsize)
    while fitted >= min_fontsize:
        temp.set_fontsize(fitted)
        bbox_disp = temp.get_window_extent(renderer=renderer)
        bbox_data = bbox_disp.transformed(ax.transData.inverted())
        if bbox_data.width <= max_width and bbox_data.height <= max_height:
            break
        fitted -= 0.5

    temp.remove()
    return max(fitted, min_fontsize)


# Helper function to draw multiline text with explicit newline handling
def draw_multiline_text(
    ax,
    x,
    y,
    text,
    fontsize,
    fontweight,
    linespacing=1.2,
    max_width=None,
    max_height=None
):
    normalized = text.replace('\\n', '\n')
    fitted_fontsize = fontsize
    if max_width is not None and max_height is not None:
        fitted_fontsize = fit_text_fontsize(
            ax,
            x,
            y,
            normalized,
            fontweight,
            fontsize,
            max_width=max_width,
            max_height=max_height,
            linespacing=linespacing
        )

    ax.text(
        x,
        y,
        normalized,
        ha='center',
        va='center',
        fontsize=fitted_fontsize,
        fontweight=fontweight,
        multialignment='center',
        linespacing=linespacing,
        zorder=25
    )

# Helper function to create a diamond shape for decision points
def draw_diamond(
    ax,
    center,
    width,
    height,
    text,
    fill_color='#E6F3FF',
    fontweight='bold',
    fontsize=None,
    fit_width_ratio=0.62,
    fit_height_ratio=0.62,
    linespacing=1.2
):
    x, y = center
    points = [
        (x, y + height/2),
        (x + width/2, y),
        (x, y - height/2),
        (x - width/2, y)
    ]
    polygon = patches.Polygon(points, closed=True, edgecolor='black', facecolor=fill_color, linewidth=2.5, zorder=10)
    ax.add_patch(polygon)
    draw_multiline_text(
        ax,
        x,
        y,
        text,
        fontsize=fontsize if fontsize is not None else DECISION_BASE_FONTSIZE,
        fontweight=fontweight,
        max_width=width * fit_width_ratio,
        max_height=height * fit_height_ratio,
        linespacing=linespacing
    )
    return {
        'n': (x, y + height/2),
        'e': (x + width/2, y),
        's': (x, y - height/2),
        'w': (x - width/2, y),
        'center': center
    }

# Helper function to create a rounded box for showing removed loops
def draw_box(ax, center, width, height, text, fill_color='white', edgecolor='black', fontweight='bold'):
    x, y = center
    corner = (x - width/2, y - height/2)
    box = patches.FancyBboxPatch(corner, width, height, boxstyle="round,pad=0.02,rounding_size=0.1", 
                                 edgecolor=edgecolor, facecolor=fill_color, linewidth=2.5, zorder=10)
    ax.add_patch(box)
    draw_multiline_text(
        ax,
        x,
        y,
        text,
        fontsize=BASE_FONTSIZE,
        fontweight=fontweight,
        max_width=width * 0.90,
        max_height=height * 0.84
    )
    return {
        'n': (x, y + height/2),
        'e': (x + width/2, y),
        's': (x, y - height/2),
        'w': (x - width/2, y),
        'center': center
    }

# Helper function to create a parallelogram for input/output data
def draw_parallelogram(ax, center, width, height, text, fill_color='white', edgecolor='black', fontweight='bold', skew=3):
    x, y = center
    points = [
        (x - width/2 + skew, y + height/2), # Top Left
        (x + width/2 + skew, y + height/2), # Top Right
        (x + width/2 - skew, y - height/2), # Bottom Right
        (x - width/2 - skew, y - height/2)  # Bottom Left
    ]
    polygon = patches.Polygon(points, closed=True, edgecolor=edgecolor, facecolor=fill_color, linewidth=2.5, zorder=10)
    ax.add_patch(polygon)
    draw_multiline_text(
        ax,
        x,
        y,
        text,
        fontsize=BASE_FONTSIZE,
        fontweight=fontweight,
        max_width=(width - 2 * abs(skew)) * 0.92,
        max_height=height * 0.84
    )
    
    return {
        'n': (x + skew, y + height/2), 
        'e': (x + width/2, y),         
        's': (x - skew, y - height/2), 
        'w': (x - width/2, y),         
        'center': center,
        'bottom_mid': (x, y - height/2),
        'top_mid': (x, y + height/2)
    }

# Helper function to draw arrows between nodes
def draw_arrow(ax, start, end, label=None, label_pos=0.5, offset_label=(0,0)):
    ax.annotate("", xy=end, xycoords='data', xytext=start, textcoords='data',
                arrowprops=dict(arrowstyle="-|>", lw=2.5, color='black'), zorder=5)
    
    if label:
        lx = start[0] + (end[0] - start[0]) * label_pos + offset_label[0]
        ly = start[1] + (end[1] - start[1]) * label_pos + offset_label[1]
        draw_multiline_text(ax, lx, ly, label, fontsize=16, fontweight='bold')

# Helper function to draw arrows with multiple segments
def draw_polyline_arrow(ax, points, label=None, label_idx=0, offset_label=(0,0)):
    for i in range(len(points)-1):
        p_start = points[i]
        p_end = points[i+1]
        if i == len(points) - 2:
            ax.annotate("", xy=p_end, xycoords='data', xytext=p_start, textcoords='data',
                        arrowprops=dict(arrowstyle="-|>", lw=2.5, color='black'), zorder=5)
        else:
             ax.plot([p_start[0], p_end[0]], [p_start[1], p_end[1]], color='black', lw=2.5, zorder=5)
    if label:
        p_s = points[label_idx]
        p_e = points[label_idx+1]
        mid_x = (p_s[0] + p_e[0])/2
        mid_y = (p_s[1] + p_e[1])/2
        draw_multiline_text(ax, mid_x + offset_label[0], mid_y + offset_label[1], label, fontsize=16, fontweight='bold')

# Setup the plot
fig, ax = plt.subplots(figsize=(20, 22)) 
ax.set_xlim(0, 178)
ax.set_ylim(5, 138)
ax.axis('off')

# Define coordinates for the nodes
y_start = 130
y_dec1 = 115
y_dec2 = 100
y_tss1 = 80
y_ctcf = 75
y_tss2 = 60
y_final = 35
y_end = 15

x_center = 75
x_left = 35
x_right = 125
x_rem_r1 = 125
x_rem_inner = 75
x_rem_outer = 165
x_rem_comb = 125

# Create the nodes of the flowchart

# Starting point: Total loops from 10 samples
n_start = draw_parallelogram(ax, (x_center, y_start), 40, 9, "58,992 loops annotated\nfrom 10 samples")

# First filter: Remove redundant loops
n_dec1 = draw_diamond(
    ax,
    (x_center, y_dec1),
    40,
    9,
    "Unique loops only",
    fontsize=26,
    fit_width_ratio=0.76,
    fit_height_ratio=0.72,
    linespacing=1.0
)
n_rem1 = draw_box(ax, (x_rem_r1, y_dec1), 28, 7, "27,219 loops\nremoved", fill_color='#FFE6E6', fontweight='normal')

# Second filter: Filter by loop length
n_dec2 = draw_diamond(
    ax,
    (x_center, y_dec2),
    40,
    9,
    "Loop length < 2 Mb",
    fontsize=26,
    fit_width_ratio=0.76,
    fit_height_ratio=0.72,
    linespacing=1.0
)
n_rem2 = draw_box(ax, (x_rem_r1, y_dec2), 26, 7, "754 loops\nremoved", fill_color='#FFE6E6', fontweight='normal')

# Parallel branches for CTCF and TSS/Promoter filtering
# CTCF branch
n_dec_ctcf = draw_diamond(
    ax,
    (x_left, y_ctcf),
    42,
    12,
    "Both anchors have\nCTCF sites (score > 6)",
    fontsize=26,
    fit_width_ratio=0.76,
    fit_height_ratio=0.72,
    linespacing=1.0
)
n_rem_ctcf = draw_box(ax, (x_rem_inner, y_ctcf), 22, 7, "5,751 loops\nremoved", fill_color='#FFE6E6', fontweight='normal')

# TSS/Promoter branch
n_dec_tss1 = draw_diamond(
    ax,
    (x_right, y_tss1),
    44,
    12,
    "TSS/Promoter assigned\nto at least one anchor",
    fontsize=26,
    fit_width_ratio=0.76,
    fit_height_ratio=0.72,
    linespacing=1.0
)
n_rem_tss1 = draw_box(ax, (x_rem_outer, y_tss1), 20, 9, "13,042\nloops\nremoved", fill_color='#FFE6E6', fontweight='normal')

n_dec_tss2 = draw_diamond(
    ax,
    (x_right, y_tss2),
    44,
    12,
    "Anchor-TSS/Promoter\ndistance < 200 kb",
    fontsize=26,
    fit_width_ratio=0.76,
    fit_height_ratio=0.72,
    linespacing=1.0
)
n_rem_tss2 = draw_box(ax, (x_rem_outer, y_tss2), 20, 9, "328\nloops\nremoved", fill_color='#FFE6E6', fontweight='normal')

# Final intersection of both branches
n_dec_final = draw_diamond(ax, (x_center, y_final), 40, 10, "Intersection of \n both criteria")
n_rem_final_combined = draw_box(
    ax,
    (x_rem_comb, y_final),
    48,
    14,
    "10,183 loops (failed CTCF)\n2,563 loops (failed TSS/Promoter:\n1,794/769) removed",
    fill_color='#FFE6E6',
    fontweight='normal'
)

# Final result
n_end = draw_parallelogram(
    ax,
    (x_center, y_end),
    50,
    13,
    "15,085 loops retained\n(CTCF-TSS: 10,125,\nCTCF-Promoter: 4,960)",
    fill_color='#E6FFCC',
    edgecolor='black'
)

# Connect the nodes with arrows and labels

# From start to first decision
draw_arrow(ax, n_start['bottom_mid'], n_dec1['n'])

# First decision: Redundant loops
draw_arrow(ax, n_dec1['e'], n_rem1['w'], label="No", offset_label=(0, 1.5))
draw_arrow(ax, n_dec1['s'], n_dec2['n'], label="Yes   (31,773 retained)", offset_label=(7.5, 0))

# Second decision: Loop length
draw_arrow(ax, n_dec2['e'], n_rem2['w'], label="No", offset_label=(0, 1.5))
y_fork_line = 90 
draw_polyline_arrow(ax, [n_dec2['s'], (x_center, y_fork_line)], label="Yes   (31,019 retained)", label_idx=0, offset_label=(7.5, 0))

# Fork to parallel branches
draw_polyline_arrow(ax, [(x_center, y_fork_line), (x_left, y_fork_line), n_dec_ctcf['n']])
draw_polyline_arrow(ax, [(x_center, y_fork_line), (x_right, y_fork_line), n_dec_tss1['n']])

# CTCF branch logic
draw_arrow(ax, n_dec_ctcf['e'], n_rem_ctcf['w'], label="No", offset_label=(0, 1.5))
y_merge = 48
draw_polyline_arrow(ax, [n_dec_ctcf['s'], (x_left, y_merge), (x_center-3.0, y_merge), (x_center-3.0, n_dec_final['n'][1])], label="Yes   (25,268 retained)", label_idx=0, offset_label=(7.5, 0)) 

# TSS/Promoter branch logic
draw_arrow(ax, n_dec_tss1['e'], n_rem_tss1['w'], label="No", offset_label=(0, 1.5))
draw_arrow(ax, n_dec_tss1['s'], n_dec_tss2['n'], label="Yes   (17,976 retained)", offset_label=(7.5, 0))
draw_arrow(ax, n_dec_tss2['e'], n_rem_tss2['w'], label="No", offset_label=(0, 1.5))
draw_polyline_arrow(ax, [n_dec_tss2['s'], (x_right, y_merge), (x_center+3.0, y_merge), (x_center+3.0, n_dec_final['n'][1])], label="Yes   (17,648 retained)\n(TSS/Promoter:11,919/5,729)", label_idx=0, offset_label=(7.5, 0))

# Final intersection logic
draw_arrow(ax, n_dec_final['e'], n_rem_final_combined['w'], label="No", offset_label=(0, 1.5))

# Final arrow to the end result
draw_arrow(ax, n_dec_final['s'], n_end['top_mid'], label="Yes", offset_label=(3.0, 0))

# Save the flowchart to a file
plt.tight_layout()
out_path_script = Path(__file__).resolve().parent / "flowchart_v11.png"
out_path_cwd = Path.cwd() / "flowchart_v11.png"
plt.savefig(out_path_script, dpi=300, bbox_inches='tight')
if out_path_cwd != out_path_script:
    plt.savefig(out_path_cwd, dpi=300, bbox_inches='tight')
print(f"Flowchart saved as {out_path_script}")
if out_path_cwd != out_path_script:
    print(f"Flowchart also saved as {out_path_cwd}")
