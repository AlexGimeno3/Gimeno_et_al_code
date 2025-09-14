import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


def create_risk_categories_table(fis_dict, base_save_dir, file_name="table_2_risk_categories.pdf"):
    tables_dir = os.path.join(base_save_dir, "Tables")
    os.makedirs(tables_dir, exist_ok=True)

    analyzed = [f for f in fis_dict.values() if f.analyzed_whole_sample == 1]
    if not analyzed:
        print("No analyzed FIS objects found.")
        return

    nd_group = [f for f in analyzed if f.ND == 1]
    non_nd_group = [f for f in analyzed if f.ND == 0]

    total_n, nd_n, non_nd_n = len(analyzed), len(nd_group), len(non_nd_group)

    sts_result = process_risk_variable(analyzed, nd_group, non_nd_group, "sts_eacts_category", [1, 2, 3, 4, 5])
    rachs_result = process_risk_variable(analyzed, nd_group, non_nd_group, "rachs1_category", [1, 2, 3, 4, 5, 6])

    create_risk_table_pdf(
        sts_result, rachs_result, total_n, nd_n, non_nd_n,
        os.path.join(tables_dir, file_name)
    )

    print(f"Table 2 saved to {os.path.join(tables_dir, file_name)}")


def process_risk_variable(all_group, nd_group, non_nd_group, attr, categories):
    def safe_get(objs):
        return [getattr(f, attr) for f in objs if getattr(f, attr) is not None]

    def count(vals, cat):
        return sum(1 for v in vals if v == cat)

    all_vals = safe_get(all_group)
    nd_vals = safe_get(nd_group)
    non_nd_vals = safe_get(non_nd_group)

    result = {"categories": categories, "all": [], "nd": [], "non_nd": []}

    for cat in categories:
        result["all"].append((count(all_vals, cat), len(all_vals)))
        result["nd"].append((count(nd_vals, cat), len(nd_vals)))
        result["non_nd"].append((count(non_nd_vals, cat), len(non_nd_vals)))

    contingency = []
    for i, cat in enumerate(categories):
        nd_c = result["nd"][i][0]
        non_c = result["non_nd"][i][0]
        if nd_c > 0 or non_c > 0:
            contingency.append([nd_c, non_c])

    try:
        chi2, p_chi, _, _ = stats.chi2_contingency(contingency)
        obs = np.array(contingency)
        total = obs.sum()
        expected = np.outer(obs.sum(1), obs.sum(0)) / total
        with np.errstate(divide="ignore", invalid="ignore"):
            lr_stat = 2 * np.nansum(obs * np.log(obs / expected))
        df = (obs.shape[0] - 1) * (obs.shape[1] - 1)
        p_lr = stats.chi2.sf(lr_stat, df)
    except Exception:
        p_chi, p_lr = None, None

    result["p_chi"] = p_chi
    result["p_lr"] = p_lr
    return result


def create_risk_table_pdf(sts, rachs, total_n, nd_n, non_nd_n, path):
    fig_height = 9
    fig, ax = plt.subplots(figsize=(10, fig_height))
    ax.axis("off")

    header = ["Category", f"Total (n={total_n})", f"ND (n={nd_n})", f"No ND (n={non_nd_n})"]
    table_data = []

    def add_section(name, result):
        all_n = result["all"][0][1]
        nd_n = result["nd"][0][1]
        non_nd_n = result["non_nd"][0][1]
        header = [
            "Category",
            f"Total (n={all_n})",
            f"ND (n={nd_n})",
            f"No ND (n={non_nd_n})"
        ]
        table_data.append([f"{name}", "", "", ""])
        table_data.append(header)
        for i, cat in enumerate(result["categories"]):
            def fmt(val):
                c, n = val
                pct = (c / n * 100) if n > 0 else 0
                return f"{c} ({pct:.1f}%)"

            table_data.append([
                str(cat),
                fmt(result["all"][i]),
                fmt(result["nd"][i]),
                fmt(result["non_nd"][i])
            ])
        p_chi = result.get("p_chi")
        p_lr = result.get("p_lr")
        if p_chi is not None and p_lr is not None:
            p_chi_str = "<0.001" if p_chi < 0.001 else f"{p_chi:.3f}"
            p_lr_str = "<0.001" if p_lr < 0.001 else f"{p_lr:.3f}"
            if p_chi < 0.05 or p_lr < 0.05:
                p_fmt = rf"$\mathbf{{{p_chi_str}/{p_lr_str}}}$"
            else:
                p_fmt = rf"${p_chi_str}/{p_lr_str}$"
        else:
            p_fmt = "N/A"
        table_data.append([r"$\chi^2$/LR test", "", "", p_fmt])

    add_section("STS-EACTS Category", sts)
    add_section("RACHS-1 Category", rachs)

    table = ax.table(
        cellText=table_data,
        colLabels=header,
        loc="center",
        cellLoc="center",
        colWidths=[0.25, 0.25, 0.25, 0.25]
    )

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.4)

    for i, row in enumerate(table_data):
        if row[0] in ["Category", "STS-EACTS Category", "RACHS-1 Category"]:
            continue
        if "Total (n=" in row[1]:
            for j in range(len(row)):
                cell = table[i + 1, j]
                cell.set_text_props(weight="bold", style="italic")

    # Style
    for i, row in enumerate(table_data):
        for j in range(len(header)):
            cell = table[i + 1, j]
            cell.visible_edges = "BRTL"
            if row[0].strip() in ["STS-EACTS Category", "RACHS-1 Category"]:
                if j == 0:
                    cell.set_text_props(weight="bold", style="italic")
                elif j == 3:
                    cell.get_text().set_text("")
                    cell.visible_edges = "BTR"
                else:
                    cell.get_text().set_text("")
                    cell.visible_edges = ""
            elif row[0] == r"$\chi^2$/LR test":
                cell.set_facecolor("0.9")

    # Title and footnote
    plt.figtext(0.05, 0.97, r"$\bf{Table\ 2:}$ Surgical risk categories", ha="left", fontsize=9.5)
    footnote = (
        "NB: χ² and LR tests compare ND vs non-ND groups. "
        "Likelihood ratio (LR) provided where expected values < 5."
    )
    plt.figtext(0.05, 0.02, footnote, ha="left", fontsize=8)

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)

def generate_risk_categories_table(stat_gui_instance):
    """
    Wrapper to call create_risk_categories_table from stat_gui.py

    Parameters:
    -----------
    stat_gui_instance : stat_gui
        An instance of the stat_gui class
    """
    base_save_dir = getattr(stat_gui_instance, 'base_save_dir', None)
    if base_save_dir is None:
        from PDF_fxns import get_base_save_dir
        base_save_dir = get_base_save_dir()

    create_risk_categories_table(stat_gui_instance.fis_dict, base_save_dir)