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
    linespacing=1.15
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
        y,
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


def draw_box_with_bullets(
    ax,
    center,
    width,
    height,
    header_text,
    bullet_items,
    fill_color='#FFE6E6',
    edgecolor='black',
    header_fontsize=13,
    bullet_fontsize=12,
    linespacing=1.25
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

    # Header text at upper center
    ax.text(
        x,
        y + height / 2 - 3.2,
        header_text,
        ha='center',
        va='top',
        fontsize=header_fontsize,
        fontweight='normal',
        color='black',
        zorder=25
    )

    # Bullet items aligned perfectly to the same left vertical line
    bullet_text = "\n".join([f"•  {item}" for item in bullet_items])
    ax.text(
        x - width * 0.40,
        y - 1.2,
        bullet_text,
        ha='left',
        va='center',
        fontsize=bullet_fontsize,
        fontweight='normal',
        linespacing=linespacing,
        color='black',
        zorder=25
    )

    return {
        'n': (x, y + height / 2),
        'e': (x + width / 2, y),
        's': (x, y - height / 2),
        'w': (x - width / 2, y),
        'center': center
    }


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
    title_fontsize=14,
    item_fontsize=12,
    linespacing=1.35
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

    # Title at upper center
    ax.text(
        x,
        y + height / 2 - 2.5,
        title,
        ha='center',
        va='top',
        fontsize=title_fontsize,
        fontweight='bold',
        color='black',
        zorder=25
    )

    # Separator horizontal line
    line_w = width * 0.78
    ax.plot([x - line_w / 2, x + line_w / 2], [y + height / 2 - 4.8, y + height / 2 - 4.8], color='#4D7C0F', lw=1.2, zorder=25)

    # Numbered list aligned to exactly the same left vertical line
    full_body = "\n".join(list_items)
    if footer_note:
        full_body += "\n" + footer_note

    ax.text(
        x - width * 0.39,
        y - 2.2,
        full_body,
        ha='left',
        va='center',
        fontsize=item_fontsize,
        fontweight='normal',
        linespacing=linespacing,
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
        max_width=width * 0.84,
        max_height=height * 0.84,
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
    fig, ax = plt.subplots(figsize=(16, 22))
    ax.set_xlim(0, 160)
    ax.set_ylim(0, 220)
    ax.axis('off')

    x_center = 75
    x_rem_right = 136
    x_left = 38
    x_right = 112

    # Y-coordinates
    y_start = 208
    y_dec1 = 188
    y_dec2 = 164
    y_dec3 = 138

    y_atac_single = 90
    y_dual_detail = 90
    y_end = 16

    # 1. Start Node
    n_start = draw_parallelogram(
        ax,
        (x_center, y_start),
        72,
        10,
        "59,000 pooled HiCCUPS loop calls\n(10 rat strains across 5K, 10K, 25K resolutions)",
        fill_color='#F0F4F8',
        edgecolor='#334E68',
        fontsize=15
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
        fontsize=15,
        fit_text=False
    )
    n_rem1 = draw_box(
        ax,
        (x_rem_right, y_dec1),
        32,
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
        fontsize=15,
        fit_text=False
    )
    n_rem2 = draw_box(
        ax,
        (x_rem_right, y_dec2),
        32,
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
        fontsize=15,
        fit_text=False
    )
    n_rem3 = draw_box(
        ax,
        (x_rem_right, y_dec3),
        36,
        11,
        "14,678 loops (47.3%)\nwithout direct\npromoter/TSS overlap",
        fill_color='#FFE6E6',
        edgecolor='black',
        fontweight='normal',
        linespacing=1.1
    )

    # 5. Branch Split: Single Promoter vs Dual Promoter
    # Left Branch: Single Direct Promoter (12,295)
    n_dec_single_atac = draw_diamond(
        ax,
        (x_left, y_atac_single),
        54,
        15,
        "Distal candidate anchor\nhas non-TSS ATAC (≥ 50 bp)?",
        fill_color='#E0F2FE',
        edgecolor='#0284C7',
        fontsize=15,
        fit_text=False
    )

    n_box_single_no_atac = draw_box(
        ax,
        (x_left - 18, y_atac_single - 24),
        34,
        12,
        "1,826 loops (5.9%)\nSingle promoter w/o\ndistal ATAC support",
        fill_color='#FFE6E6',
        edgecolor='black',
        fontweight='normal',
        linespacing=1.1
    )

    n_box_single_putative = draw_box(
        ax,
        (x_left + 18, y_atac_single - 24),
        34,
        12,
        "10,469 loops (33.7%)\nPutative Regulatory Loops\n(Open-chromatin P-E)",
        fill_color='#FFE6E6',
        edgecolor='black',
        fontweight='normal',
        linespacing=1.1
    )

    # Right Branch: Dual Promoter (4,048)
    n_box_dual = draw_box_with_bullets(
        ax,
        (x_right, y_dual_detail),
        54,
        22,
        header_text="4,048 Both-Anchor Promoter Loops (13.0%)",
        bullet_items=[
            "P–E: 1,094 (27.0%)",
            "P–(P/E): 629 (15.5%)",
            "(P/E)–(P/E): 1,184 (29.2%)",
            "No distal E: 1,141 (28.2%)"
        ],
        fill_color='#FFE6E6',
        edgecolor='black',
        header_fontsize=13,
        bullet_fontsize=12,
        linespacing=1.25
    )

    # 6. Final Comprehensive Master Summary
    n_final_master = draw_parallelogram_summary(
        ax,
        (x_center, y_end + 12),
        138,
        22,
        title="Final Categorization of 31,021 Pooled Chromatin Loops (< 2 Mb)",
        list_items=[
            "1. Putative Regulatory Loops (promoter/TSS + distal non-TSS ATAC): 10,469 (33.7%)",
            "2. Both-anchor Promoter Loops (promoter/TSS direct overlap at both anchors): 4,048 (13.0%)",
            "3. Single-anchor Promoter Loops without distal ATAC support: 1,826 (5.9%)",
            "4. Structural Loops without direct promoter/TSS overlap: 14,678 (47.3%)"
        ],
        footer_note="* Independent Structural Annotation: 29,980 loops (96.6%) have predicted CTCF motifs at both anchors",
        fill_color='#ECFCCB',
        edgecolor='#4D7C0F',
        title_fontsize=14,
        item_fontsize=12,
        linespacing=1.35
    )

    # Connecting Arrows
    draw_arrow(ax, n_start['bottom_mid'], n_dec1['n'])
    draw_arrow(ax, n_dec1['e'], n_rem1['w'], label="No", offset_label=(0, 1.5))
    draw_arrow(ax, n_dec1['s'], n_dec2['n'], label="Yes (31,778 retained)", offset_label=(0, 0))

    draw_arrow(ax, n_dec2['e'], n_rem2['w'], label="No", offset_label=(0, 1.5))
    draw_arrow(ax, n_dec2['s'], n_dec3['n'], label="Yes (31,021 distinct loops)", offset_label=(0, 0))

    draw_arrow(ax, n_dec3['e'], n_rem3['w'], label="None", offset_label=(0, 1.5))

    # Fork from Dec3 to Single vs Dual
    y_split_pt = 122
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
        [(x_center, y_split_pt), (x_right, y_split_pt), n_box_dual['n']],
        label="Direct overlap at\nboth anchors (4,048 loops)",
        label_idx=0,
        offset_label=(0, 2.5),
        label_fontsize=11
    )

    # ATAC decision arrows (Left Branch) - labels placed at vertical segment midpoints
    draw_polyline_arrow(
        ax,
        [n_dec_single_atac['w'], (x_left - 18, n_dec_single_atac['w'][1]), n_box_single_no_atac['n']],
        label="No  (< 50 bp)",
        label_idx=1,
        offset_label=(0, 0),
        label_fontsize=10
    )
    draw_polyline_arrow(
        ax,
        [n_dec_single_atac['e'], (x_left + 18, n_dec_single_atac['e'][1]), n_box_single_putative['n']],
        label="Yes  (≥ 50 bp)",
        label_idx=1,
        offset_label=(0, 0),
        label_fontsize=10
    )

    # Arrows to Final Summary Node
    draw_arrow(ax, n_box_single_no_atac['s'], (x_left - 18, n_final_master['top_mid'][1]))
    draw_arrow(ax, n_box_single_putative['s'], (x_left + 18, n_final_master['top_mid'][1]))
    draw_arrow(ax, n_box_dual['s'], (x_right, n_final_master['top_mid'][1]))

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
    out_base = results_dir / "flowchart_revision"
    generate_revision_flowchart(out_base)
