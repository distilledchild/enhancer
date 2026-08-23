import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path

# Base typography configuration
BASE_FONTSIZE = 14
MIN_FONTSIZE = 9
DECISION_BASE_FONTSIZE = 16


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


def draw_rich_line_centered(ax, x, y, text, fontsize, fontweight, color='black', scale_badge=1.25):
    """Draw a single line centered horizontally at (x, y), scaling ❶, ❷, ❸, ❹, ❺ by scale_badge with exact metrics."""
    CIRCLED_SYMBOLS = {'❶', '❷', '❸', '❹', '❺', '①', '②', '③', '④', '⑤'}
    has_symbol = any(s in text for s in CIRCLED_SYMBOLS)
    if not has_symbol:
        ax.text(x, y, text, ha='center', va='center', fontsize=fontsize, fontweight=fontweight, color=color, zorder=25)
        return

    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    # Find the symbol
    symbol = [s for s in CIRCLED_SYMBOLS if s in text][0]
    idx = text.index(symbol)
    part_before = text[:idx]
    part_symbol = symbol
    part_after = text[idx + 1:]

    # Measure exact widths using sentinel characters
    def measure_token(txt, fsize, fweight):
        if not txt:
            return 0
        t_flank = ax.text(0, 0, f"X{txt}X", fontsize=fsize, fontweight=fweight, alpha=0.0)
        t_base = ax.text(0, 0, "XX", fontsize=fsize, fontweight=fweight, alpha=0.0)
        w = (
            t_flank.get_window_extent(renderer=renderer).transformed(ax.transData.inverted()).width
            - t_base.get_window_extent(renderer=renderer).transformed(ax.transData.inverted()).width
        )
        t_flank.remove()
        t_base.remove()
        return w

    w_sp = measure_token(" ", fontsize, fontweight)

    has_space_after = part_after.startswith(' ')
    stripped_after = part_after.lstrip(' ')

    w_before = measure_token(part_before, fontsize, fontweight) if part_before else 0
    t_symbol = ax.text(0, 0, part_symbol, fontsize=fontsize * scale_badge, fontweight='bold', alpha=0.0)
    w_symbol = t_symbol.get_window_extent(renderer=renderer).transformed(ax.transData.inverted()).width
    t_symbol.remove()
    w_after = measure_token(stripped_after, fontsize, fontweight) if stripped_after else 0

    total_w = w_before + w_symbol + (w_sp if has_space_after else 0) + w_after
    start_x = x - total_w / 2

    cur_x = start_x
    if part_before:
        ax.text(cur_x, y, part_before, ha='left', va='center', fontsize=fontsize, fontweight=fontweight, color=color, zorder=25)
        cur_x += w_before

    ax.text(cur_x, y, part_symbol, ha='left', va='center', fontsize=fontsize * scale_badge, fontweight='bold', color=color, zorder=25)
    cur_x += w_symbol

    if has_space_after:
        cur_x += w_sp

    if stripped_after:
        ax.text(cur_x, y, stripped_after, ha='left', va='center', fontsize=fontsize, fontweight=fontweight, color=color, zorder=25)


def draw_multiline_text(
    ax,
    x,
    y,
    text,
    fontsize,
    fontweight,
    linespacing=1.2,
    max_width=None,
    max_height=None,
    color='black'
):
    normalized = text.replace('\\n', '\n')
    lines = normalized.split('\n')
    n_lines = len(lines)

    CIRCLED_SYMBOLS = {'❶', '❷', '❸', '❹', '❺', '①', '②', '③', '④', '⑤'}
    has_any_symbol = any(any(s in l for s in CIRCLED_SYMBOLS) for l in lines)

    if not has_any_symbol:
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
            color=color,
            multialignment='center',
            linespacing=linespacing,
            zorder=25
        )
        return

    # Draw line by line to support enlarged badge
    line_step = fontsize * 0.28 * linespacing
    total_h = (n_lines - 1) * line_step
    top_y = y + total_h / 2

    for i, line in enumerate(lines):
        line_y = top_y - i * line_step
        draw_rich_line_centered(ax, x, line_y, line, fontsize, fontweight, color=color, scale_badge=1.25)


def draw_diamond(
    ax,
    center,
    width,
    height,
    text,
    fill_color='#E0F2FE',
    edgecolor='#0284C7',
    fontweight='bold',
    fontsize=None,
    fit_text=False,
    fit_width_ratio=0.75,
    fit_height_ratio=0.75,
    linespacing=1.15,
    y_offset=0.0
):
    x, y = center
    points = [
        (x, y + height / 2),
        (x + width / 2, y),
        (x, y - height / 2),
        (x - width / 2, y)
    ]
    polygon = patches.Polygon(points, closed=True, edgecolor=edgecolor, facecolor=fill_color, linewidth=2.0, zorder=10)
    ax.add_patch(polygon)
    draw_multiline_text(
        ax,
        x,
        y + y_offset,
        text,
        fontsize=fontsize if fontsize is not None else DECISION_BASE_FONTSIZE,
        fontweight=fontweight,
        max_width=width * fit_width_ratio if fit_text else None,
        max_height=height * fit_height_ratio if fit_text else None,
        linespacing=linespacing
    )
    return {
        'n': (x, y + height / 2),
        'e': (x + width / 2, y),
        's': (x, y - height / 2),
        'w': (x - width / 2, y),
        'center': center
    }


def draw_box(
    ax,
    center,
    width,
    height,
    text,
    fill_color='white',
    edgecolor='#555555',
    fontweight='bold',
    fontsize=None,
    fit_text=True,
    linespacing=1.15
):
    x, y = center
    corner = (x - width / 2, y - height / 2)
    box = patches.FancyBboxPatch(
        corner,
        width,
        height,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        edgecolor=edgecolor,
        facecolor=fill_color,
        linewidth=2.0,
        zorder=10
    )
    ax.add_patch(box)
    draw_multiline_text(
        ax,
        x,
        y,
        text,
        fontsize=fontsize if fontsize is not None else BASE_FONTSIZE,
        fontweight=fontweight,
        max_width=width * 0.90 if fit_text else None,
        max_height=height * 0.90 if fit_text else None,
        linespacing=linespacing
    )
    return {
        'n': (x, y + height / 2),
        'e': (x + width / 2, y),
        's': (x, y - height / 2),
        'w': (x - width / 2, y),
        'center': center
    }


def draw_rich_line_left(ax, x, y, text, fontsize, fontweight, color='black', scale_badge=1.20):
    """Draw a left-aligned line with scaled ❶, ❷, ❸, ❹, with exact font metrics and spacing."""
    CIRCLED_SYMBOLS = {'❶', '❷', '❸', '❹', '①', '②', '③', '④'}
    has_symbol = any(s in text for s in CIRCLED_SYMBOLS)
    if not has_symbol:
        ax.text(x, y, text, ha='left', va='center', fontsize=fontsize, fontweight=fontweight, color=color, zorder=25)
        return

    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    import re
    tokens = re.split(r'([❶❷❸❹①②③④])', text)
    cur_x = x

    def measure_token(txt, fsize, fweight):
        if not txt:
            return 0
        t_flank = ax.text(0, 0, f"X{txt}X", fontsize=fsize, fontweight=fweight, alpha=0.0)
        t_base = ax.text(0, 0, "XX", fontsize=fsize, fontweight=fweight, alpha=0.0)
        w = (
            t_flank.get_window_extent(renderer=renderer).transformed(ax.transData.inverted()).width
            - t_base.get_window_extent(renderer=renderer).transformed(ax.transData.inverted()).width
        )
        t_flank.remove()
        t_base.remove()
        return w

    for token in tokens:
        if not token:
            continue
        if token in CIRCLED_SYMBOLS:
            t = ax.text(cur_x, y, token, ha='left', va='center', fontsize=fontsize * scale_badge, fontweight='bold', color=color, zorder=25)
            w = t.get_window_extent(renderer=renderer).transformed(ax.transData.inverted()).width
            cur_x += w
        else:
            l_spaces = len(token) - len(token.lstrip(' '))
            r_spaces = len(token) - len(token.rstrip(' '))
            core_text = token.strip(' ')

            if l_spaces > 0:
                cur_x += measure_token(' ' * l_spaces, fontsize, fontweight)

            if core_text:
                ax.text(cur_x, y, core_text, ha='left', va='center', fontsize=fontsize, fontweight=fontweight, color=color, zorder=25)
                cur_x += measure_token(core_text, fontsize, fontweight)

            if r_spaces > 0:
                cur_x += measure_token(' ' * r_spaces, fontsize, fontweight)


def draw_parallelogram_summary(
    ax,
    center,
    width,
    height,
    title,
    list_items,
    footer_note=None,
    fill_color='#ECFCCB',
    edgecolor='#4D7C0F',
    title_fontsize=14.5,
    item_fontsize=14.5,
    linespacing=1.25
):
    x, y = center
    shift = width * 0.05
    points = [
        (x - width / 2 + shift, y + height / 2),
        (x + width / 2 + shift, y + height / 2),
        (x + width / 2 - shift, y - height / 2),
        (x - width / 2 - shift, y - height / 2)
    ]
    polygon = patches.Polygon(points, closed=True, edgecolor=edgecolor, facecolor=fill_color, linewidth=2.5, zorder=10)
    ax.add_patch(polygon)

    # Title at upper center
    ax.text(
        x,
        y + height / 2 - 2.4,
        title,
        ha='center',
        va='top',
        fontsize=title_fontsize,
        fontweight='bold',
        color='black',
        zorder=25
    )

    # Separator horizontal line
    line_w = width * 0.86
    sep_y = y + height / 2 - 4.6
    ax.plot([x - line_w / 2, x + line_w / 2], [sep_y, sep_y], color='#4D7C0F', lw=1.2, zorder=25)

    # Numbered list aligned to left vertical line
    all_lines = list_items[:]
    if footer_note:
        all_lines.append(footer_note)

    n_lines = len(all_lines)
    # Vertically center the text block between sep_y and bottom of parallelogram (y - height/2)
    bottom_y = y - height / 2
    avail_height = sep_y - bottom_y
    line_step = 2.70
    block_height = (n_lines - 1) * line_step
    top_y = sep_y - (avail_height - block_height) / 2

    for i, line_str in enumerate(all_lines):
        cur_y = top_y - i * line_step
        # Bold the main category headers (1., 2., 3.)
        is_header = line_str.strip().startswith(('1.', '2.', '3.'))
        line_weight = 'bold' if is_header else 'normal'
        ax.text(
            x - width * 0.42,
            cur_y,
            line_str,
            ha='left',
            va='center',
            fontsize=item_fontsize,
            fontweight=line_weight,
            color='black',
            zorder=25
        )

    return {
        'top_mid': (x, y + height / 2),
        'bottom_mid': (x, y - height / 2),
        'center': center
    }


def draw_parallelogram(
    ax,
    center,
    width,
    height,
    text,
    fill_color='#EBF7E3',
    edgecolor='#2E7D32',
    fontweight='bold',
    fontsize=None,
    fit_text=False,
    linespacing=1.2
):
    x, y = center
    shift = width * 0.06
    points = [
        (x - width / 2 + shift, y + height / 2),
        (x + width / 2 + shift, y + height / 2),
        (x + width / 2 - shift, y - height / 2),
        (x - width / 2 - shift, y - height / 2)
    ]
    polygon = patches.Polygon(points, closed=True, edgecolor=edgecolor, facecolor=fill_color, linewidth=2.5, zorder=10)
    ax.add_patch(polygon)
    draw_multiline_text(
        ax,
        x,
        y,
        text,
        fontsize=fontsize if fontsize is not None else BASE_FONTSIZE + 2,
        fontweight=fontweight,
        max_width=width * 0.84 if fit_text else None,
        max_height=height * 0.84 if fit_text else None,
        linespacing=linespacing
    )
    return {
        'top_mid': (x, y + height / 2),
        'bottom_mid': (x, y - height / 2),
        'center': center
    }


def draw_arrow(ax, start, end, label=None, offset_label=(0, 0), label_fontsize=11, label_color='black', color='black', lw=2.0):
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            lw=lw,
            mutation_scale=16,
            shrinkA=0,
            shrinkB=0
        ),
        zorder=5
    )
    if label:
        x_mid = (start[0] + end[0]) / 2 + offset_label[0]
        y_mid = (start[1] + end[1]) / 2 + offset_label[1]
        ax.text(
            x_mid,
            y_mid,
            label,
            ha='center',
            va='center',
            fontsize=label_fontsize,
            fontweight='bold',
            color=label_color,
            bbox=dict(boxstyle="square,pad=0.15", facecolor="white", edgecolor="none", alpha=0.9),
            zorder=20
        )


def draw_polyline_arrow(ax, points, label=None, label_idx=0, offset_label=(0, 0), label_fontsize=11, color='black', lw=2.0):
    for i in range(len(points) - 1):
        is_last = (i == len(points) - 2)
        p1 = points[i]
        p2 = points[i + 1]
        if is_last:
            ax.annotate(
                "",
                xy=p2,
                xytext=p1,
                arrowprops=dict(
                    arrowstyle="-|>",
                    color=color,
                    lw=lw,
                    mutation_scale=16,
                    shrinkA=0,
                    shrinkB=0
                ),
                zorder=5
            )
        else:
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, lw=lw, zorder=5)

        if label and i == label_idx:
            x_mid = (p1[0] + p2[0]) / 2 + offset_label[0]
            y_mid = (p1[1] + p2[1]) / 2 + offset_label[1]
            ax.text(
                x_mid,
                y_mid,
                label,
                ha='center',
                va='center',
                fontsize=label_fontsize,
                fontweight='bold',
                bbox=dict(boxstyle="square,pad=0.15", facecolor="white", edgecolor="none", alpha=0.9),
                zorder=20
            )


def generate_revision_flowchart(output_file):
    fig, ax = plt.subplots(figsize=(17, 22))
    ax.set_xlim(-15, 178)
    ax.set_ylim(0, 220)
    ax.axis('off')

    x_center = 74
    x_rem_right = 152  # Shifted right to create ample space for vertical alignment
    x_left = 37
    x_right = 107

    # Y-coordinates
    y_start = 208
    y_dec1 = 188
    y_dec2 = 164
    y_dec3 = 138

    y_atac_split = 92
    y_leaf_boxes = 58
    y_end = 16

    # 1. Start Node
    n_start = draw_parallelogram(
        ax,
        (x_center, y_start),
        76,
        11.5,
        "59,000 pooled HiCCUPS loop calls\n(10 rat strains across 5K, 10K, 25K resolutions)",
        fill_color='#F0F4F8',
        edgecolor='#334E68',
        fontsize=14.5,
        fit_text=False,
        linespacing=1.2
    )

    # 2. Decision 1: Distinct Coordinates
    n_dec1 = draw_diamond(
        ax,
        (x_center, y_dec1),
        48,
        13,
        "Collapse to unique\nloop coordinates",
        fill_color='#E0F2FE',
        edgecolor='#0284C7',
        fontsize=14.5,
        fit_text=False
    )
    n_rem1 = draw_box(
        ax,
        (x_rem_right, y_dec1),
        34,
        9,
        "27,222 replicate calls\ncollapsed",
        fill_color='#FFE6E6',
        edgecolor='black',
        fontweight='normal'
    )

    # 3. Decision 2: Distance < 2 Mb (Reference Benchmark)
    n_dec2 = draw_diamond(
        ax,
        (x_center, y_dec2),
        48,
        13,
        "Loop distance < 2 Mb\n(Intrachromosomal)",
        fill_color='#E0F2FE',
        edgecolor='#0284C7',
        fontsize=14.5,
        fit_text=False
    )
    n_rem2 = draw_box(
        ax,
        (x_rem_right, y_dec2),
        34,
        9,
        "757 loops (≥ 2 Mb)\nexcluded",
        fill_color='#FFE6E6',
        edgecolor='black',
        fontweight='normal'
    )

    # 4. Decision 3: Direct Promoter/TSS Overlap
    n_dec3 = draw_diamond(
        ax,
        (x_center, y_dec3),
        56,
        14,
        "Direct Promoter/TSS overlap\n(TSS ± 1 kb at ≥ 1 anchor)",
        fill_color='#E0F2FE',
        edgecolor='#0284C7',
        fontsize=14.5,
        fit_text=False
    )
    # Box 14,678 with black border (Set to 13.0pt as requested)
    n_rem3 = draw_box(
        ax,
        (x_rem_right, y_dec3),
        36,
        13,
        "14,678 loops (47.3%)\n❺ Loops without\ndirect promoter/TSS overlap",
        fill_color='#FFE6E6',
        edgecolor='black',
        fontsize=13.0,
        fontweight='normal',
        linespacing=1.12
    )

    # Consistent font size for bottom 4 leaf boxes (12.0pt to fit perfectly)
    LEAF_FONTSIZE = 12.0

    # 5. Branch Split: Single Promoter vs Dual Promoter
    # Left Branch: Single Direct Promoter (12,295 loops)
    n_dec_single_atac = draw_diamond(
        ax,
        (x_left, y_atac_split),
        52,
        14,
        "Distal candidate anchor\nhas non-TSS ATAC (≥ 50 bp)?",
        fill_color='#E0F2FE',
        edgecolor='#0284C7',
        fontsize=14.5,
        fit_text=False,
        y_offset=-0.2
    )

    # Left = Yes (Green), Right = No (Orange)
    # Box 1: Single Putative (❶)
    n_box_single_putative = draw_box(
        ax,
        (19.5, y_leaf_boxes),
        32,
        15,
        "10,469 loops (33.7%)\n❶ Single-promoter\nPutative Regulatory",
        fill_color='#DCFCE7',
        edgecolor='#16A34A',
        fontsize=LEAF_FONTSIZE,
        fontweight='normal',
        fit_text=False,
        linespacing=1.12
    )

    # Box 2: Single No ATAC (❷)
    n_box_single_no_atac = draw_box(
        ax,
        (53.0, y_leaf_boxes),
        32,
        15,
        "1,826 loops (5.9%)\n❷ Single-promoter w/o\ndistal ATAC support",
        fill_color='#FFF7ED',
        edgecolor='#C2410C',
        fontsize=LEAF_FONTSIZE,
        fontweight='normal',
        fit_text=False,
        linespacing=1.12
    )

    # Right Branch: Both-Anchor Promoter (4,048 loops)
    n_dec_dual_atac = draw_diamond(
        ax,
        (x_right, y_atac_split),
        54,
        14,
        "Residual non-TSS ATAC (≥ 50 bp)\nat ≥ 1 anchor?",
        fill_color='#E0F2FE',
        edgecolor='#0284C7',
        fontsize=14.5,
        fit_text=False,
        y_offset=-0.75
    )

    # Box 3: Dual Putative (❸) - width=39.0, center=90.0
    x_box3_center = 90.0
    n_box_dual_putative = draw_box(
        ax,
        (x_box3_center, y_leaf_boxes),
        39.0,
        15,
        "2,907 loops (9.4%)\n❸ Dual-promoter with\nEnhancer ATAC Support\n(P–E:1,094 | P–P/E:629 | P/E–P/E:1,184)",
        fill_color='#DCFCE7',
        edgecolor='#16A34A',
        fontsize=LEAF_FONTSIZE,
        fontweight='normal',
        fit_text=False,
        linespacing=1.08
    )

    # Box 4: Dual No ATAC (❹) - center=128.0, width=34.0
    x_box4_center = 128.0
    n_box_dual_no_atac = draw_box(
        ax,
        (x_box4_center, y_leaf_boxes),
        34.0,
        15,
        "1,141 loops (3.7%)\n❹ Pure Promoter–\nPromoter Contacts\n(No distal enhancer ATAC)",
        fill_color='#FFF7ED',
        edgecolor='#C2410C',
        fontsize=LEAF_FONTSIZE,
        fontweight='normal',
        fit_text=False,
        linespacing=1.12
    )

    # 6. Final Comprehensive Master Summary (Expanded horizontally to catch vertical line at x=152)
    n_final_master = draw_parallelogram_summary(
        ax,
        (x_center + 4, y_end + 11.5),
        156,
        25,
        title="Final Categorization of 31,021 Pooled Chromatin Loops (< 2 Mb)",
        list_items=[
            "1. Putative regulatory loops with distal open-chromatin support: 13,376 loops (43.1%)",
            "   • ❶ Single-promoter with distal ATAC (10,469) + ❸ Dual-promoter with enhancer ATAC (2,907)",
            "2. Promoter-associated loops without distal ATAC support: 2,967 loops (9.6%)",
            "   • ❷ Single-promoter w/o ATAC (1,826) + ❹ Pure promoter–promoter contacts (1,141)",
            "3. Loops without direct promoter/TSS overlap: 14,678 loops (47.3%)",
            "   • ❺ No direct promoter/TSS overlap (14,678; 47.3%); 96.6% CTCF-supported structural contacts"
        ],
        # footer_note="* Independent Structural Annotation: 29,980 loops (96.6%) have predicted CTCF motifs at both anchors",
        fill_color='#ECFCCB',
        edgecolor='#4D7C0F',
        title_fontsize=14.5,
        item_fontsize=14.5,
        linespacing=1.30
    )

    # Connecting Arrows
    draw_arrow(ax, n_start['bottom_mid'], n_dec1['n'])
    draw_arrow(ax, n_dec1['e'], n_rem1['w'], label="No", offset_label=(0, 1.5))
    draw_arrow(ax, n_dec1['s'], n_dec2['n'], label="Yes (31,778 retained)", offset_label=(0, 0))

    draw_arrow(ax, n_dec2['e'], n_rem2['w'], label="No", offset_label=(0, 1.5))
    draw_arrow(ax, n_dec2['s'], n_dec3['n'], label="Yes (31,021 distinct loops)", offset_label=(0, 0))

    draw_arrow(ax, n_dec3['e'], n_rem3['w'], label="None", offset_label=(0, 1.5))

    # Fork from Dec3 to Single vs Dual
    y_split_pt = 120
    draw_polyline_arrow(
        ax,
        [n_dec3['s'], (x_center, y_split_pt), (x_left, y_split_pt), n_dec_single_atac['n']],
        label="Direct overlap at\nsingle anchor (12,295 loops)",
        label_idx=1,
        offset_label=(0, 2.5),
        label_fontsize=11
    )
    draw_polyline_arrow(
        ax,
        [n_dec3['s'], (x_center, y_split_pt), (x_right, y_split_pt), n_dec_dual_atac['n']],
        label="Direct overlap at\nboth anchors (4,048 loops)",
        label_idx=1,
        offset_label=(0, 2.5),
        label_fontsize=11
    )

    # ATAC decision arrows (Left Branch: Left=Yes, Right=No)
    draw_polyline_arrow(
        ax,
        [n_dec_single_atac['w'], (x_left - 17.5, n_dec_single_atac['w'][1]), n_box_single_putative['n']],
        label="Yes  (≥ 50 bp)",
        label_idx=1,
        offset_label=(0, 0),
        label_fontsize=10
    )
    draw_polyline_arrow(
        ax,
        [n_dec_single_atac['e'], (53.0, n_dec_single_atac['e'][1]), n_box_single_no_atac['n']],
        label="No  (< 50 bp)",
        label_idx=1,
        offset_label=(0, 0),
        label_fontsize=10
    )

    # ATAC decision arrows (Right Branch: Left=Yes, Right=No)
    draw_polyline_arrow(
        ax,
        [n_dec_dual_atac['w'], (x_box3_center, n_dec_dual_atac['w'][1]), n_box_dual_putative['n']],
        label="Yes  (≥ 50 bp)",
        label_idx=1,
        offset_label=(0, 0),
        label_fontsize=10
    )
    draw_polyline_arrow(
        ax,
        [n_dec_dual_atac['e'], (x_box4_center, n_dec_dual_atac['e'][1]), n_box_dual_no_atac['n']],
        label="No  (< 50 bp)",
        label_idx=1,
        offset_label=(0, 0),
        label_fontsize=10
    )

    # Arrows to Final Summary Node
    # Category 1 arrows (Green)
    draw_arrow(ax, n_box_single_putative['s'], (x_left - 17.5, n_final_master['top_mid'][1]), color='#16A34A', lw=2.2)
    draw_arrow(ax, n_box_dual_putative['s'], (x_box3_center, n_final_master['top_mid'][1]), color='#16A34A', lw=2.2)

    # Category 2 arrows (Orange/Dark)
    draw_arrow(ax, n_box_single_no_atac['s'], (x_left + 16.5, n_final_master['top_mid'][1]), color='#C2410C', lw=2.0)
    draw_arrow(ax, n_box_dual_no_atac['s'], (x_box4_center, n_final_master['top_mid'][1]), color='#C2410C', lw=2.0)

    # Category 3: Perfectly straight vertical black line from Box 14,678 down into the parallelogram!
    draw_arrow(
        ax,
        n_rem3['s'],
        (x_rem_right, n_final_master['top_mid'][1]),
        color='black',
        lw=2.0
    )

    plt.tight_layout()
    output_png = Path(output_file).with_suffix(".png")
    output_pdf = Path(output_file).with_suffix(".pdf")
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.savefig(output_pdf, bbox_inches='tight')
    plt.close()
    print(f"Revision flowchart successfully saved to:\n  - PNG: {output_png}\n  - PDF: {output_pdf}")


if __name__ == "__main__":
    results_dir = Path(__file__).resolve().parent / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    out_base = results_dir / "revision_figure7_flowchart"
    generate_revision_flowchart(out_base)
