import matplotlib.pyplot as plt
import matplotlib.patches as patches

# This script generates a flowchart illustrating the filtering process of chromatin loops.
# It uses matplotlib to draw shapes and text, representing each step of the analysis.

# Helper function to create a diamond shape for decision points
def draw_diamond(ax, center, width, height, text, fill_color='#E6F3FF', fontweight='bold'):
    x, y = center
    points = [
        (x, y + height/2),
        (x + width/2, y),
        (x, y - height/2),
        (x - width/2, y)
    ]
    polygon = patches.Polygon(points, closed=True, edgecolor='black', facecolor=fill_color, linewidth=1.5, zorder=10)
    ax.add_patch(polygon)
    ax.text(x, y, text, ha='center', va='center', fontsize=9, wrap=True, zorder=20, fontweight=fontweight)
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
                                 edgecolor=edgecolor, facecolor=fill_color, linewidth=1.5, zorder=10)
    ax.add_patch(box)
    ax.text(x, y, text, ha='center', va='center', fontsize=9, wrap=True, zorder=20, fontweight=fontweight)
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
    polygon = patches.Polygon(points, closed=True, edgecolor=edgecolor, facecolor=fill_color, linewidth=1.5, zorder=10)
    ax.add_patch(polygon)
    ax.text(x, y, text, ha='center', va='center', fontsize=9, wrap=True, zorder=20, fontweight=fontweight)
    
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
                arrowprops=dict(arrowstyle="-|.", lw=1.5, color='black'), zorder=5)
    
    if label:
        lx = start[0] + (end[0] - start[0]) * label_pos + offset_label[0]
        ly = start[1] + (end[1] - start[1]) * label_pos + offset_label[1]
        ax.text(lx, ly, label, ha='center', va='center', fontsize=9, zorder=30, fontweight='bold',
                    bbox=dict(facecolor='white', edgecolor='none', pad=2))

# Helper function to draw arrows with multiple segments
def draw_polyline_arrow(ax, points, label=None, label_idx=0, offset_label=(0,0)):
    for i in range(len(points)-1):
        p_start = points[i]
        p_end = points[i+1]
        if i == len(points) - 2:
            ax.annotate("", xy=p_end, xycoords='data', xytext=p_start, textcoords='data',
                        arrowprops=dict(arrowstyle="-|.", lw=1.5, color='black'), zorder=5)
        else:
             ax.plot([p_start[0], p_end[0]], [p_start[1], p_end[1]], color='black', lw=1.5, zorder=5)
    if label:
        p_s = points[label_idx]
        p_e = points[label_idx+1]
        mid_x = (p_s[0] + p_e[0])/2
        mid_y = (p_s[1] + p_e[1])/2
        ax.text(mid_x + offset_label[0], mid_y + offset_label[1], label, ha='center', va='center', fontsize=9, zorder=30, fontweight='bold',
                bbox=dict(facecolor='white', edgecolor='none', pad=2))

# Setup the plot
fig, ax = plt.subplots(figsize=(12, 14)) 
ax.set_xlim(0, 110)
ax.set_ylim(20, 100)
ax.axis('off')

# Define coordinates for the nodes
y_start = 96
y_dec1 = 88
y_dec2 = 80
y_tss1 = 68
y_ctcf = 60
y_tss2 = 56
y_final = 42
y_end = 30

x_center = 55
x_left = 28
x_right = 82
x_rem_r1 = 90
x_rem_inner = 55
x_rem_outer = 105
x_rem_comb = 90

# Create the nodes of the flowchart

# Starting point: Total loops from 10 samples
n_start = draw_parallelogram(ax, (x_center, y_start), 26, 5, "58,992 loops annotated\nfrom 10 samples")

# First filter: Remove redundant loops
n_dec1 = draw_diamond(ax, (x_center, y_dec1), 26, 5, "Unique loops only")
n_rem1 = draw_box(ax, (x_rem_r1, y_dec1), 18, 4, "27,219 loops\nremoved", fill_color='#FFE6E6', fontweight='normal')

# Second filter: Filter by loop length
n_dec2 = draw_diamond(ax, (x_center, y_dec2), 26, 5, "Loop length < 2 Mb")
n_rem2 = draw_box(ax, (x_rem_r1, y_dec2), 16, 4, "754 loops\nremoved", fill_color='#FFE6E6', fontweight='normal')

# Parallel branches for CTCF and TSS/Promoter filtering
# CTCF branch
n_dec_ctcf = draw_diamond(ax, (x_left, y_ctcf), 28, 7, "Both anchors have\nCTCF sites (score > 6)")
n_rem_ctcf = draw_box(ax, (x_rem_inner, y_ctcf), 14, 4, "5,751 loops\nremoved", fill_color='#FFE6E6', fontweight='normal')

# TSS/Promoter branch
n_dec_tss1 = draw_diamond(ax, (x_right, y_tss1), 30, 7, "TSS/Promoter assigned to\n at least one anchor")
n_rem_tss1 = draw_box(ax, (x_rem_outer, y_tss1), 12, 5, "13,042\nloops\nremoved", fill_color='#FFE6E6', fontweight='normal')

n_dec_tss2 = draw_diamond(ax, (x_right, y_tss2), 30, 7, "Anchor-TSS/Promoter distance\n< 3rd Quartile")
n_rem_tss2 = draw_box(ax, (x_rem_outer, y_tss2), 12, 5, "1,864\nloops\nremoved", fill_color='#FFE6E6', fontweight='normal')

# Final intersection of both branches
n_dec_final = draw_diamond(ax, (x_center, y_final), 26, 6, "Intersection of both criteria")
n_rem_final_combined = draw_box(ax, (x_rem_comb, y_final), 28, 6, "11,445 loops (failed CTCF)\n2,286 loops (failed TSS/Promoter)\nremoved", fill_color='#FFE6E6', fontweight='normal')

# Final result
n_end = draw_parallelogram(ax, (x_center, y_end), 24, 5, "13,823 loops\nretained", fill_color='#E6FFCC', edgecolor='black')

# Connect the nodes with arrows and labels

# From start to first decision
draw_arrow(ax, n_start['bottom_mid'], n_dec1['n'])

# First decision: Redundant loops
draw_arrow(ax, n_dec1['e'], n_rem1['w'], label="No", offset_label=(0, 0.8))
draw_arrow(ax, n_dec1['s'], n_dec2['n'], label="Yes (31,773 retained)", offset_label=(4.5, 0))

# Second decision: Loop length
draw_arrow(ax, n_dec2['e'], n_rem2['w'], label="No", offset_label=(0, 0.8))
y_fork_line = 74 
draw_polyline_arrow(ax, [n_dec2['s'], (x_center, y_fork_line)], label="Yes (31,019 retained)", label_idx=0, offset_label=(4.5, 0))

# Fork to parallel branches
draw_polyline_arrow(ax, [(x_center, y_fork_line), (x_left, y_fork_line), n_dec_ctcf['n']])
draw_polyline_arrow(ax, [(x_center, y_fork_line), (x_right, y_fork_line), n_dec_tss1['n']])

# CTCF branch logic
draw_arrow(ax, n_dec_ctcf['e'], n_rem_ctcf['w'], label="No", offset_label=(0, 0.8))
y_merge = 48
draw_polyline_arrow(ax, [n_dec_ctcf['s'], (x_left, y_merge), (x_center-2, y_merge), (x_center-2, n_dec_final['n'][1])], label="Yes (25,268 retained)", label_idx=0, offset_label=(4.5, 0)) 

# TSS/Promoter branch logic
draw_arrow(ax, n_dec_tss1['e'], n_rem_tss1['w'], label="No", offset_label=(0, 0.8))
draw_arrow(ax, n_dec_tss1['s'], n_dec_tss2['n'], label="Yes (17,977 retained)", offset_label=(4.5, 0))
draw_arrow(ax, n_dec_tss2['e'], n_rem_tss2['w'], label="No", offset_label=(0, 0.8))
draw_polyline_arrow(ax, [n_dec_tss2['s'], (x_right, y_merge), (x_center+2, y_merge), (x_center+2, n_dec_final['n'][1])], label="Yes (16,113 retained)", label_idx=0, offset_label=(4.5, 0))

# Final intersection logic
draw_arrow(ax, n_dec_final['e'], n_rem_final_combined['w'], label="No", offset_label=(0, 0.8))

# Final arrow to the end result
draw_arrow(ax, n_dec_final['s'], n_end['top_mid'], label="Yes", offset_label=(2.0, 0))

# Save the flowchart to a file
plt.tight_layout()
plt.savefig('flowchart_v12.png', dpi=300, bbox_inches='tight')
print("Flowchart saved as flowchart_v12.png")
