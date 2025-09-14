import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import re

def create_neurodevelopmental_table(fis_dict, base_save_dir, file_name="table_4_neurodevelopmental_outcomes.pdf"):
    """
    Create and save a neurodevelopmental outcomes table comparing ND and non-ND groups.
    
    Parameters:
    -----------
    fis_dict : dict
        Dictionary of FIS objects
    base_save_dir : str
        Base directory to save the table
    file_name : str
        Name of the output PDF file
    """
    #Step 0: Create directory to save table
    tables_dir = os.path.join(base_save_dir, "Tables")
    if not os.path.exists(tables_dir):
        os.makedirs(tables_dir)
    
    #Step 1: Filter for analyzed cases only. Quality check that filters correctly reflect the analyzed patients.
    analyzed_fis = [fis_obj for fis_obj in fis_dict.values() if fis_obj.analyzed_whole_sample == 1]
    filter_ok_fis = [fis_obj for fis_obj in fis_dict.values() if fis_obj.filter == 1]
    if not len(filter_ok_fis)==len(analyzed_fis):
        raise ValueError(f"{len(analyzed_fis)} patients were analyzed, but {len(filter_ok_fis)} were let through via filtering.")
    
    if not analyzed_fis:
        print("No analyzed FIS objects found. Cannot create neurodevelopmental outcomes table.")
        return
    
    # Split into groups by ND status
    all_group = analyzed_fis
    nd_group = [fis_obj for fis_obj in analyzed_fis if (fis_obj.ND == 1 and not fis_obj.ND is None)]
    non_nd_group = [fis_obj for fis_obj in analyzed_fis if (fis_obj.ND == 0 and not fis_obj.ND is None)]
    if not (len(all_group)== len(nd_group)+len(non_nd_group)):
        raise ValueError(f"The analyzed group had {len(all_group)} patients, but there were {len(nd_group)} in the ND group and {len(non_nd_group)} in the non-ND group.")
    
    # Total counts
    total_n = len(all_group) #Total n is number of patients with the analyzed flag as 1
    nd_n = len(nd_group) #nd_n is number of analyzed patients in the ND group
    non_nd_n = len(non_nd_group) #non_nd_n is number of analyzed patients in the non-ND group
    
    #print(f"Total analyzed cases: {total_n}")
    #print(f"ND group: {nd_n}")
    #print(f"Non-ND group: {non_nd_n}")

    # Get counts for 'ND present?'
    nd_present_count = len(nd_group)
    
    # Define table structure
    table_data = []
    
    # Add ND present count
    table_data.append({
        "name": "ND present?",
        "is_header": True,
        "is_section_header": True,
        "all": f"{nd_present_count}/{total_n} ({(nd_present_count/total_n*100):.1f}%)",
        "nd": "",
        "non_nd": "",
        "p_value": "",
        "p_value_float": None,
        "test_symbol": ""
    })
    
    # Add Bayley's section header
    table_data.append({
        "name": "Bayley's characteristics",
        "is_header": False,
        "is_section_header": True,
        "all": "",
        "nd": "",
        "non_nd": "",
        "p_value": "",
        "p_value_float": None,
        "test_symbol": ""
    })
    
    # Get Bayley test counts and filter objects by test type
    bayley_all = [obj for obj in all_group if obj.test_used == 2]
    bayley_nd = [obj for obj in nd_group if obj.test_used == 2]
    bayley_non_nd = [obj for obj in non_nd_group if obj.test_used == 2]
    
    bayley_n = len(bayley_all)
    bayley_nd_n = len(bayley_nd)
    bayley_non_nd_n = len(bayley_non_nd)
    
    # Add Bayley's variables
    bayley_vars = [
        {"name": "Cognitive composite score", "attr": "bayley_cognitive_composite"},
        {"name": "Language composite score", "attr": "bayley_language_composite"},
        {"name": "Motor composite score", "attr": "bayley_motor_composite"},
        {"name": "Days of life at test", "attr": "dol_for_test_used"}
    ]
    
    for var in bayley_vars:
        result = process_continuous_var(
            var["name"], 
            [getattr(obj, var["attr"]) for obj in bayley_all if getattr(obj, var["attr"]) is not None],
            [getattr(obj, var["attr"]) for obj in bayley_nd if getattr(obj, var["attr"]) is not None],
            [getattr(obj, var["attr"]) for obj in bayley_non_nd if getattr(obj, var["attr"]) is not None],
            "auto"
        )
        
        # Update the sample counts in the section headers
        result["all"] = result["all"].replace("sample", f"sample (n = {bayley_n})")
        result["nd"] = result["nd"].replace("ND", f"ND (n = {bayley_nd_n})")
        result["non_nd"] = result["non_nd"].replace("No ND", f"No ND (n = {bayley_non_nd_n})")
        
        table_data.append(result)
    
    # Add Vineland section header
    table_data.append({
        "name": "Vineland characteristics",
        "is_header": False,
        "is_section_header": True,
        "all": "",
        "nd": "",
        "non_nd": "",
        "p_value": "",
        "p_value_float": None,
        "test_symbol": ""
    })
    
    # Get Vineland test counts and filter objects by test type
    vineland_all = [obj for obj in all_group if obj.test_used == 1]
    vineland_nd = [obj for obj in nd_group if obj.test_used == 1]
    vineland_non_nd = [obj for obj in non_nd_group if obj.test_used == 1]
    
    vineland_n = len(vineland_all)
    vineland_nd_n = len(vineland_nd)
    vineland_non_nd_n = len(vineland_non_nd)
    
    # Add Vineland variables
    vineland_vars = [
        {"name": "Communication total score", "attr": "vl_communication"},
        {"name": "Activities of daily living total score", "attr": "vl_daily_activities"},
        {"name": "Socialization total score", "attr": "vl_socialization"},
        {"name": "Development total score", "attr": "vl_psychomotor_development"},
        {"name": "Adaptive behavior composite score", "attr": "vl_abc"},
        {"name": "Days of life at test", "attr": "dol_for_test_used"}
    ]
    
    for var in vineland_vars:
        result = process_continuous_var(
            var["name"], 
            [getattr(obj, var["attr"]) for obj in vineland_all if getattr(obj, var["attr"]) is not None],
            [getattr(obj, var["attr"]) for obj in vineland_nd if getattr(obj, var["attr"]) is not None],
            [getattr(obj, var["attr"]) for obj in vineland_non_nd if getattr(obj, var["attr"]) is not None],
            "auto"
        )
        
        # Update the sample counts in the section headers
        result["all"] = result["all"].replace("sample", f"sample (n = {vineland_n})")
        result["nd"] = result["nd"].replace("ND", f"ND (n = {vineland_nd_n})")
        result["non_nd"] = result["non_nd"].replace("No ND", f"No ND (n = {vineland_non_nd_n})")
        
        table_data.append(result)
    
    # Create the PDF table
    create_table_pdf(table_data, total_n, nd_n, non_nd_n, bayley_n, bayley_nd_n, bayley_non_nd_n, vineland_n, vineland_nd_n, vineland_non_nd_n, os.path.join(tables_dir, file_name))
    
    print(f"Neurodevelopmental outcomes table saved to {os.path.join(tables_dir, file_name)}")
    return table_data

def process_continuous_var(name, all_data, nd_data, non_nd_data, test_type):
    """Process a continuous variable and decide appropriate test"""
    result = {
        "name": name, 
        "is_header": False,
        "is_section_header": False
    }
    p_value = None
    p_value_float = None
    test_symbol = ""

    # If no data for either group, return empty result
    if not nd_data or not non_nd_data:
        return {
            "name": name,
            "is_header": False,
            "is_section_header": False,
            "all": "Insufficient data",
            "nd": "Insufficient data",
            "non_nd": "Insufficient data",
            "p_value": "-",
            "p_value_float": None,
            "test_symbol": ""
        }

    # Decide which test to use if auto
    if test_type == "auto":
        # Ensure data is not empty before shapiro test
        if len(nd_data) >= 3 and len(non_nd_data) >= 3 and len(np.unique(nd_data)) > 1 and len(np.unique(non_nd_data)) > 1:
            try:
                _, p_nd = stats.shapiro(nd_data)
                _, p_non_nd = stats.shapiro(non_nd_data)
                if p_nd > 0.05 and p_non_nd > 0.05:
                    test_type = "ttest"
                else:
                    test_type = "mwu"
            except ValueError:  # Handle potential errors in shapiro test
                test_type = "mwu"
        else:
            test_type = "mwu"  # Default to non-parametric for small or constant samples

    # Get common stats only if data exists
    mean_all, median_all, q1_all, q3_all = np.nan, np.nan, np.nan, np.nan
    if all_data:
        all_data = [x for x in all_data if x is not None and not np.isnan(x)]
        if all_data:
            mean_all = np.mean(all_data)
            median_all = np.median(all_data)
            q1_all = np.percentile(all_data, 25)
            q3_all = np.percentile(all_data, 75)

    mean_nd, median_nd, q1_nd, q3_nd = np.nan, np.nan, np.nan, np.nan
    if nd_data:
        nd_data = [x for x in nd_data if x is not None and not np.isnan(x)]
        if nd_data:
            mean_nd = np.mean(nd_data)
            median_nd = np.median(nd_data)
            q1_nd = np.percentile(nd_data, 25)
            q3_nd = np.percentile(nd_data, 75)

    mean_non_nd, median_non_nd, q1_non_nd, q3_non_nd = np.nan, np.nan, np.nan, np.nan
    if non_nd_data:
        non_nd_data = [x for x in non_nd_data if x is not None and not np.isnan(x)]
        if non_nd_data:
            mean_non_nd = np.mean(non_nd_data)
            median_non_nd = np.median(non_nd_data)
            q1_non_nd = np.percentile(non_nd_data, 25)
            q3_non_nd = np.percentile(non_nd_data, 75)

    # Run the appropriate test
    valid_nd_data = [x for x in nd_data if x is not None and not np.isnan(x)]
    valid_non_nd_data = [x for x in non_nd_data if x is not None and not np.isnan(x)]
    
    if test_type == "ttest" and len(valid_nd_data) > 1 and len(valid_non_nd_data) > 1:
        try:
            _, p_var = stats.levene(valid_nd_data, valid_non_nd_data)
            equal_var = p_var >= 0.05
            t_stat, p_value = stats.ttest_ind(valid_nd_data, valid_non_nd_data, equal_var=equal_var)
            p_value_float = p_value  # Store raw float
            test_symbol = "*"

            sem_all = np.std(valid_nd_data + valid_non_nd_data, ddof=1) / np.sqrt(len(valid_nd_data) + len(valid_non_nd_data)) if valid_nd_data + valid_non_nd_data else 0
            sem_nd = np.std(valid_nd_data, ddof=1) / np.sqrt(len(valid_nd_data)) if valid_nd_data else 0
            sem_non_nd = np.std(valid_non_nd_data, ddof=1) / np.sqrt(len(valid_non_nd_data)) if valid_non_nd_data else 0

            result["all"] = f"{mean_all:.1f} ({sem_all:.1f})"
            result["nd"] = f"{mean_nd:.1f} ({sem_nd:.1f})"
            result["non_nd"] = f"{mean_non_nd:.1f} ({sem_non_nd:.1f})"
        except (ValueError, ZeroDivisionError):
            p_value = np.nan
            p_value_float = np.nan
            test_symbol = "*"  # Still indicate t-test was attempted
            result["all"] = "N/A"
            result["nd"] = "N/A"
            result["non_nd"] = "N/A"
    else:  # MWU
        if valid_nd_data and valid_non_nd_data:
            try:
                u_stat, p_value = stats.mannwhitneyu(valid_nd_data, valid_non_nd_data, alternative='two-sided')
                p_value_float = p_value
            except ValueError:
                p_value = np.nan
                p_value_float = np.nan
        else:
            p_value = np.nan
            p_value_float = np.nan

        test_symbol = "†"
        
        # Format based on variable name for special cases
        if "Days of life" in name:
            # For Days of life, show median [min, max]
            all_min, all_max = np.min(valid_nd_data + valid_non_nd_data), np.max(valid_nd_data + valid_non_nd_data) if valid_nd_data + valid_non_nd_data else (np.nan, np.nan)
            nd_min, nd_max = np.min(valid_nd_data), np.max(valid_nd_data) if valid_nd_data else (np.nan, np.nan)
            non_nd_min, non_nd_max = np.min(valid_non_nd_data), np.max(valid_non_nd_data) if valid_non_nd_data else (np.nan, np.nan)
            
            result["all"] = f"{median_all:.0f} [{all_min:.0f}, {all_max:.0f}]" if not np.isnan(median_all) else "N/A"
            result["nd"] = f"{median_nd:.0f} [{nd_min:.0f}, {nd_max:.0f}]" if not np.isnan(median_nd) else "N/A"
            result["non_nd"] = f"{median_non_nd:.0f} [{non_nd_min:.0f}, {non_nd_max:.0f}]" if not np.isnan(median_non_nd) else "N/A"
        else:
            # For other variables, show median [Q1, Q3]
            result["all"] = f"{median_all:.1f} [{q1_all:.1f}, {q3_all:.1f}]" if not np.isnan(median_all) else "N/A"
            result["nd"] = f"{median_nd:.1f} [{q1_nd:.1f}, {q3_nd:.1f}]" if not np.isnan(median_nd) else "N/A"
            result["non_nd"] = f"{median_non_nd:.1f} [{q1_non_nd:.1f}, {q3_non_nd:.1f}]" if not np.isnan(median_non_nd) else "N/A"

    # Format p-value for table display
    if p_value is not None and not np.isnan(p_value):
        if p_value < 0.001:
            result["p_value"] = "< 0.001"
        else:
            result["p_value"] = f"{p_value:.3f}"
    else:
        result["p_value"] = "-"
    
    result["p_value_float"] = p_value_float
    result["test_symbol"] = test_symbol
    
    return result

def create_table_pdf(table_data, total_n, nd_n, non_nd_n, bayley_n, bayley_nd_n, bayley_non_nd_n, vineland_n, vineland_nd_n, vineland_non_nd_n, output_path):
    """Creates and saves a formatted neurodevelopmental outcomes table PDF."""
    # --- Data Preparation ---
    table_rows = []
    p_value_metadata = []

    # First row: ND present?
    first_row = table_data[0]
    table_rows.append([first_row["name"], first_row["all"], "", "", ""])
    p_value_metadata.append(None)

    # Build table content
    for row_data in table_data[1:]:
        if row_data.get("is_section_header", False):
            name = row_data["name"]
            if "Bayley" in name:
                table_rows.append([
                    name,
                    f"Total sample (n = {bayley_n})",
                    f"ND (n = {bayley_nd_n})",
                    f"No ND (n = {bayley_non_nd_n})",
                    ""
                ])
            elif "Vineland" in name:
                table_rows.append([
                    name,
                    f"Total sample (n = {vineland_n})",
                    f"ND (n = {vineland_nd_n})",
                    f"No ND (n = {vineland_non_nd_n})",
                    ""
                ])
            else:
                table_rows.append([name, "Total sample", "ND", "No ND", ""])
            p_value_metadata.append(None)
        else:
            # Clean out sample size annotations from text
            all_val = re.sub(r'\s*\(n\s*=\s*\d+\)', '', str(row_data["all"]))
            nd_val = re.sub(r'\s*\(n\s*=\s*\d+\)', '', str(row_data["nd"]))
            non_nd_val = re.sub(r'\s*\(n\s*=\s*\d+\)', '', str(row_data["non_nd"]))
            p_val_str = row_data.get("p_value", "")
            p_val_float = row_data.get("p_value_float", None)
            test_symbol = row_data.get("test_symbol", "")

            table_rows.append([row_data["name"], all_val, nd_val, non_nd_val, p_val_str])
            p_value_metadata.append({'float': p_val_float, 'symbol': test_symbol, 'p_value': p_val_str})

    # --- Plot Settings ---
    base_fontsize = 8.5
    row_height = 0.28
    header_height = 0.5
    title_space = 0.3
    footer_space = 1.5
    fig_height = len(table_rows) * row_height + header_height + title_space + footer_space

    fig, ax = plt.subplots(figsize=(12, fig_height))
    ax.axis('off')

    # Layout margins
    left_margin = 0.05
    right_margin = 0.95
    bottom_margin = 0.15
    top_margin = 0.95
    title_y_offset = 0.92
    footnote_y_start = 0.08
    footnote_spacing = 0.035

    plt.subplots_adjust(left=left_margin, right=right_margin, bottom=bottom_margin, top=top_margin)

    # Create table
    table = ax.table(
        cellText=table_rows,
        loc='center',
        cellLoc='left',
        colWidths=[0.35, 0.23, 0.20, 0.20, 0.07],
    )

    # Format cells
    for i in range(len(table_rows)):
        for j in range(5):
            try:
                cell = table[(i, j)]
            except KeyError:
                continue

            # Determine if this is a section header row
            is_section = table_rows[i][0] in ["Bayley's characteristics", "Vineland characteristics"]

            if is_section:
                cell.set_text_props(weight='bold', style='italic')
                if j == 0:
                    cell.visible_edges = 'LR'
                elif j == 4:
                    cell.visible_edges = 'LR'
                else:
                    cell.visible_edges = 'LR'  # Remove inner vertical lines
            else:
                cell.visible_edges = 'BTRL'

            cell.PAD = 0.02
            cell._loc = 'left'
            text_obj = cell.get_text()
            text_obj.set_horizontalalignment('left')
            text_obj.set_verticalalignment('center')

            # Bold significant p-values
            if j == 4 and i > 0 and p_value_metadata[i] is not None:
                p_val_float = p_value_metadata[i].get('float', None)
                if p_val_float is not None and not np.isnan(p_val_float) and p_val_float < 0.05:
                    text_obj.set_weight('bold')

    table.auto_set_font_size(False)
    table.set_fontsize(base_fontsize)
    table.scale(1, 1.5)

    # --- Title ---
    plt.figtext(left_margin, title_y_offset,
                r"$\bf{Table\ 4:}$ Neurodevelopmental outcomes data",
                ha='left', va='bottom',
                fontsize=base_fontsize)

    # --- Footnotes ---
    footnote_text = (
        "NB: Bayley's Scales of Infant and Toddler Development (3rd Ed.) and Vineland Adaptive Behavior Scales (2nd Ed.) were used; for patients who had both, the Bayley scale was used. While some patients received both tests, only the test relevant to the determination of ND was included for each patient."
    )
    test_symbols_text = (
        r"$^*$Student t$-$test used" + "\n" +
        r"$^\dagger$Mann$-$Whitney U$-$test used"
    )

    plt.figtext(left_margin, footnote_y_start, footnote_text,
                ha='left', va='top', fontsize=base_fontsize, wrap=True)
    plt.figtext(left_margin, footnote_y_start - footnote_spacing, test_symbols_text,
                ha='left', va='top', fontsize=base_fontsize)

    # Adjust bottom spacing
    lowest_footnote_pos = footnote_y_start - footnote_spacing - 0.06
    plt.subplots_adjust(bottom=max(0.01, lowest_footnote_pos))

    # Save figure
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close(fig)

    print(f"Neurodevelopmental outcomes table saved to {output_path}")

def generate_neurodevelopmental_table(stat_gui_instance):
    """
    Function to be called from stat_gui.py to generate neurodevelopmental outcomes table
    
    Parameters:
    -----------
    stat_gui_instance : stat_gui
        An instance of the stat_gui class
    """
    # Get the base save directory from the instance
    base_save_dir = getattr(stat_gui_instance, 'base_save_dir', None)
    
    if base_save_dir is None:
        # If no base_save_dir, prompt for one
        from PDF_fxns import get_base_save_dir
        base_save_dir = get_base_save_dir()
    
    # Create the demographics table
    create_neurodevelopmental_table(stat_gui_instance.fis_dict, base_save_dir)
    
    print("Neurodevelopmental outcomes table generated successfully.")

if __name__ == "__main__":
    # This can be used for testing
    try:
        from stat_gui import stat_gui
        my_gui = stat_gui()
        generate_neurodevelopmental_table(my_gui)
    except ImportError:
        print("Running in standalone mode - import a fis_dict to test")