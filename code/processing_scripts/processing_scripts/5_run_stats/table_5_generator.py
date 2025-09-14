import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def create_artefact_table(fis_dict, base_save_dir, file_name="table_5_artefact_comparison.pdf"):
    tables_dir = os.path.join(base_save_dir, "Tables")
    os.makedirs(tables_dir, exist_ok=True)

    analyzed = [f for f in fis_dict.values() if f.analyzed_whole_sample == 1]
    if not analyzed:
        print("No analyzed FIS objects found.")
        return

    nd_group = [f for f in analyzed if f.ND == 1]
    non_nd_group = [f for f in analyzed if f.ND == 0]

    total_n, nd_n, non_nd_n = len(analyzed), len(nd_group), len(non_nd_group)

    def get_values(attr):
        all_vals = [getattr(f, attr, None) for f in analyzed if getattr(f, attr, None) is not None]
        nd_vals = [getattr(f, attr, None) for f in nd_group if getattr(f, attr, None) is not None]
        non_nd_vals = [getattr(f, attr, None) for f in non_nd_group if getattr(f, attr, None) is not None]
        return all_vals, nd_vals, non_nd_vals

    def process_continuous(name, all_vals, nd_vals, non_nd_vals):
        # Skip if either group has no data
        if not nd_vals or not non_nd_vals:
            return [name, "Insufficient", "Insufficient", "Insufficient", "-"]

        # Decide which test to use
        try:
            if len(nd_vals) >= 3 and len(non_nd_vals) >= 3:
                _, p1 = stats.shapiro(nd_vals)
                _, p2 = stats.shapiro(non_nd_vals)
                normal = p1 > 0.05 and p2 > 0.05
            else:
                normal = False
        except:
            normal = False

        if normal:
            _, p_var = stats.levene(nd_vals, non_nd_vals)
            equal_var = p_var > 0.05
            stat, p = stats.ttest_ind(nd_vals, non_nd_vals, equal_var=equal_var)
            fmt = lambda x: f"{np.mean(x):.2f} (±{np.std(x, ddof=1):.2f})"
        else:
            stat, p = stats.mannwhitneyu(nd_vals, non_nd_vals, alternative='two-sided')
            fmt = lambda x: f"{np.median(x):.2f} [{np.percentile(x,25):.2f}, {np.percentile(x,75):.2f}]"

        p_str = f"{p:.3f}" if p >= 0.001 else "<0.001"
        return [name, fmt(all_vals), fmt(nd_vals), fmt(non_nd_vals), p_str]

    # Run comparisons
    rows = []
    rows.append(process_continuous("Proportion artefact", *get_values("proportion_artefact")))
    rows.append(process_continuous("Number of artefacts", *get_values("num_artefacts")))
    if any(hasattr(f, "time_artefact") for f in analyzed):
        rows.append(process_continuous("Time in artefact (hr)", *get_values("time_artefact")))

    # PDF output
    fig_height = 2 + 0.6 * len(rows)
    fig, ax = plt.subplots(figsize=(10, fig_height))
    ax.axis("off")

    columns = ["Variable", f"All (n={total_n})", f"ND (n={nd_n})", f"No ND (n={non_nd_n})", "p-value"]
    table = ax.table(
        cellText=rows,
        colLabels=columns,
        loc='center',
        cellLoc='left'
    )

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.4)

    # Format header row
    for j in range(len(columns)):
        cell = table[0, j]
        cell.set_text_props(weight='bold', style='italic', ha='left')
        cell.visible_edges = 'BRTL'

    # Format body rows
    for i, row in enumerate(rows):
        for j, val in enumerate(row):
            cell = table[i+1, j]
            cell.visible_edges = 'BRTL'
            text_obj = cell.get_text()
            text_obj.set_horizontalalignment('left')
            text_obj.set_verticalalignment('center')

            if j == 4:  # p-value column
                try:
                    p_val = float(val) if "<" not in val else 0.0001
                    if p_val < 0.05:
                        text_obj.set_weight('bold')
                except:
                    pass

    plt.figtext(0.05, 0.95, r"$\bf{Table\ 5:}$ Artefact characteristics by ND group", ha="left", fontsize=10)
    plt.figtext(0.05, 0.02, "NB: Values are mean ± SD if normal, or median [Q1, Q3] if non-normal. Mann-Whitney U or t-tests used.", ha="left", fontsize=8)

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    plt.savefig(os.path.join(tables_dir, file_name), dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Table 5 saved to {os.path.join(tables_dir, file_name)}")

def generate_artefact_table(stat_gui_instance):
    base_save_dir = getattr(stat_gui_instance, 'base_save_dir', None)
    if base_save_dir is None:
        from PDF_fxns import get_base_save_dir
        base_save_dir = get_base_save_dir()

    create_artefact_table(stat_gui_instance.fis_dict, base_save_dir)