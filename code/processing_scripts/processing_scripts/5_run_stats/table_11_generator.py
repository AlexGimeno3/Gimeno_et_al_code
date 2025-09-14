import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

def create_biomarker_table(fis_dict, base_save_dir, file_name="table_11_biomarkers.pdf"):
    tables_dir = os.path.join(base_save_dir, "Tables")
    os.makedirs(tables_dir, exist_ok=True)

    analyzed = [f for f in fis_dict.values() if f.analyzed_whole_sample == 1]
    if not analyzed:
        print("No analyzed FIS objects found.")
        return

    nd_group = [f for f in analyzed if f.ND == 1]
    non_nd_group = [f for f in analyzed if f.ND == 0]

    def extract(attr):
        return [getattr(f, attr) for f in analyzed if getattr(f, attr) is not None], \
               [getattr(f, attr) for f in nd_group if getattr(f, attr) is not None], \
               [getattr(f, attr) for f in non_nd_group if getattr(f, attr) is not None]

    def format_entry(name, all_vals, nd_vals, non_nd_vals):
        # Count entries
        counts = (len(all_vals), len(nd_vals), len(non_nd_vals))
        label = f"{name} (n = {counts[0]}, {counts[1]}, {counts[2]})"

        # Skip row if no data in either group
        if not nd_vals or not non_nd_vals:
            return [label, "Insufficient", "Insufficient", "Insufficient", "-"]

        # Decide normality
        normal = False
        try:
            if len(nd_vals) >= 3 and len(non_nd_vals) >= 3:
                _, p1 = stats.shapiro(nd_vals)
                _, p2 = stats.shapiro(non_nd_vals)
                normal = p1 > 0.05 and p2 > 0.05
        except:
            pass

        if normal:
            fmt = lambda x: f"{np.mean(x):.1f} ({np.std(x, ddof=1):.1f})"
            _, p_val = stats.ttest_ind(nd_vals, non_nd_vals, equal_var=True)
        else:
            fmt = lambda x: f"{np.median(x):.1f} [{np.percentile(x, 25):.1f}, {np.percentile(x, 75):.1f}]"
            _, p_val = stats.mannwhitneyu(nd_vals, non_nd_vals, alternative='two-sided')

        p_fmt = f"{p_val:.2f}" if p_val >= 0.01 else "<0.01"

        return [label, fmt(all_vals), fmt(nd_vals), fmt(non_nd_vals), p_fmt]

    rows = []
    rows.append(format_entry("Preop enolase", *extract("preop_enolase")))
    rows.append(format_entry("Preop S100B", *extract("preop_S100B")))
    rows.append(format_entry("Postop enolase", *extract("postop_enolase")))
    rows.append(format_entry("Postop S100B", *extract("postop_S100B")))

    fig_height = 2 + 0.6 * len(rows)
    fig, ax = plt.subplots(figsize=(9, fig_height))
    ax.axis("off")

    columns = [r"$\it{Biomarker}$", r"$\it{Total\ sample}$", r"$\it{ND}$", r"$\it{No\ ND}$", r"$\it{p\text{-}value}$"]
    table = ax.table(
        cellText=rows,
        colLabels=columns,
        loc='center',
        cellLoc='left'
    )

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.4)

    for i in range(len(rows)+1):  # +1 for header
        for j in range(len(columns)):
            cell = table[i, j]
            cell.visible_edges = 'BRTL'
            text_obj = cell.get_text()
            text_obj.set_horizontalalignment('left')
            text_obj.set_verticalalignment('center')

            # Header styling
            if i == 0:
                text_obj.set_weight('bold')
                text_obj.set_style('italic')

            # First column styling
            if j == 0:
                text_obj.set_weight('bold')

            # Bold significant p-values
            if j == 4 and i > 0:
                try:
                    val = rows[i-1][j]
                    p_val = float(val) if "<" not in val else 0.009
                    if p_val < 0.05:
                        text_obj.set_weight('bold')
                except:
                    pass

    # Manually set column widths: first wider
    for (row, col), cell in table.get_celld().items():
        if col == 0:
            cell.set_width(0.42)
        else:
            cell.set_width(0.14)

    plt.figtext(0.05, 0.95, r"$\bf{Table\ 11:}$ Biomarkers data", ha="left", fontsize=10)
    plt.figtext(0.05, 0.02,
                "NB: Values are mean (SD) if normal, or median [Q1, Q3] if non-normal. Mann-Whitney U or t-tests used.",
                ha="left", fontsize=8)

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    plt.savefig(os.path.join(tables_dir, file_name), dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Table 11 saved to {os.path.join(tables_dir, file_name)}")

def generate_biomarker_table(stat_gui_instance):
    base_save_dir = getattr(stat_gui_instance, 'base_save_dir', None)
    if base_save_dir is None:
        from PDF_fxns import get_base_save_dir
        base_save_dir = get_base_save_dir()

    create_biomarker_table(stat_gui_instance.fis_dict, base_save_dir)
