import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def create_demographics_table(fis_dict, base_save_dir, file_name="table_1_demographics.pdf", case = None):
    """
    Create and save a demographics table comparing ND and non-ND groups.
    This specific implementation focuses on the exact variables in the example table.
    
    Parameters:
    -----------
    fis_dict : dict
        Dictionary of FIS objects
    base_save_dir : str
        Base directory to save the table
    file_name : str
        Name of the output PDF file
    """
    # Create tables directory if it doesn't exist
    tables_dir = os.path.join(base_save_dir, "Tables")
    if not os.path.exists(tables_dir):
        os.makedirs(tables_dir)
    
    # Filter for analyzed cases only
    if case == "whole_sample":
        analyzed_fis = [fis_obj for fis_obj in fis_dict.values() if fis_obj.analyzed_whole_sample == 1]
    elif case == "seizure_sample":
        analyzed_fis = [fis_obj for fis_obj in fis_dict.values() if fis_obj.analyzed_seizure_sample == 1]
    elif case == "across_seizures":
        analyzed_fis = [fis_obj for fis_obj in fis_dict.values() if fis_obj.analyzed_across_seizures == 1]
    
    if not analyzed_fis:
        print("No analyzed FIS objects found. Cannot create demographics table.")
        return
    
    # Split into ND and non-ND groups
    nd_group = [fis_obj for fis_obj in analyzed_fis if fis_obj.ND == 1]
    non_nd_group = [fis_obj for fis_obj in analyzed_fis if fis_obj.ND == 0]
    
    # Print sample sizes
    total_n = len(analyzed_fis)
    nd_n = len(nd_group)
    non_nd_n = len(non_nd_group)
    
    #print(f"Total analyzed cases: {total_n}")
    #print(f"ND group: {nd_n}")
    #print(f"Non-ND group: {non_nd_n}")
    
    # Define variables to include in the table based on the provided example
    table_structure = [
        {
            "section": "Demographic Variables",
            "variables": [
                {
                    "name": "Sex (number male, n = *)",
                    "attr": "sex",
                    "type": "categorical",
                    "format": lambda x: sum(1 for obj in x if obj.sex == 2),  # Count males (sex==2)
                    "test": "chi2"
                },
                {
                    "name": "Days of life at EEG (n = *)",
                    "attr": "pma_eeg",
                    "type": "continuous",
                    "format": lambda x: [obj.pma_eeg for obj in x if obj.pma_eeg is not None],
                    "test": "auto"
                },
                {
                    "name": "Gestational age at EEG (days, n = *)",
                    "attr": "ga_eeg",
                    "type": "continuous",
                    "format": lambda x: [obj.ga_eeg for obj in x if obj.ga_eeg is not None],
                    "test": "auto"
                }
            ]
        },
        {
            "section": "Clinical variables",
            "variables": [
                {
                    "name": "Aristotle score (n = *)",
                    "attr": "aristotle_score",
                    "type": "continuous",
                    "format": lambda x: [obj.aristotle_score for obj in x if obj.aristotle_score is not None],
                    "test": "auto"
                }
            ]
        },
        {
            "section": "Surgery Details",
            "variables": [
                {
                    "name": "CPB used? (n = *)",
                    "attr": "cpb_used",
                    "type": "categorical",
                    "format": lambda x: sum(1 for obj in x if obj.cpb_used == 1),  # Count CPB used
                    "test": "chi2"
                },
                {
                    "name": "CPB time (mins, n = *)",
                    "attr": "cpb_time",
                    "type": "continuous",
                    "format": lambda x: [obj.cpb_time for obj in x if obj.cpb_time is not None],
                    "test": "auto"
                },
                {
                    "name": "Clamp used? (n = *)",
                    "attr": "cpb_clamp_used",
                    "type": "categorical",
                    "format": lambda x: sum(1 for obj in x if obj.cpb_clamp_used == 1),  # Count clamp used
                    "test": "chi2"
                },
                {
                    "name": "Time under clamp (mins, n = *)",
                    "attr": "cpb_clamp_time",
                    "type": "continuous",
                    "format": lambda x: [obj.cpb_clamp_time for obj in x if obj.cpb_clamp_time is not None],
                    "test": "auto"
                },
                {
                    "name": "Deep hypothermic circulatory arrest (DHCA) used? (n = *)",
                    "attr": "circulatory_arrest_used",
                    "type": "categorical",
                    "format": lambda x: sum(1 for obj in x if obj.circulatory_arrest_used == 1),  # Count DHCA used
                    "test": "chi2"
                },
                {
                    "name": "Time under DHCA (n = *)",
                    "attr": "circulatory_arrest_time",
                    "type": "continuous",
                    "format": lambda x: [obj.circulatory_arrest_time for obj in x if obj.circulatory_arrest_time is not None],
                    "test": "auto"
                },
                {
                    "name": "Minimum body temperature (°C, n = *)",
                    "attr": "cooling_temperature",
                    "type": "continuous",
                    "format": lambda x: [obj.cooling_temperature for obj in x if obj.cooling_temperature is not None],
                    "test": "auto"
                },
                {
                    "name": "Aortic clamp used? (n = *)",
                    "attr": "clamp_used_coarct_repair",
                    "type": "categorical",
                    "format": lambda x: sum(1 for obj in x if obj.clamp_used_coarct_repair == 1),  # Count aortic clamp used
                    "test": "chi2"
                },
                {
                    "name": "Total time under aortic clamp (n = *)",
                    "attr": "clamp_time_coarct_repair",
                    "type": "continuous",
                    "format": lambda x: [obj.clamp_time_coarct_repair for obj in x if obj.clamp_time_coarct_repair is not None],
                    "test": "auto"
                },
                {
                    "name": "Ultrafiltration used? (n = *)",
                    "attr": "uf_used",
                    "type": "categorical",
                    "format": lambda x: sum(1 for obj in x if obj.uf_used == 1),  # Count UF used
                    "test": "chi2"
                },
                {
                    "name": "Volume ultrafiltrate (mL, n = *)",
                    "attr": "uf_volume",
                    "type": "continuous",
                    "format": lambda x: [obj.uf_volume for obj in x if obj.uf_volume is not None],
                    "test": "auto"
                }
            ]
        }
    ]
    
    # Process each variable and calculate statistics
    table_data = []
    
    for section in table_structure:
        # Add section header
        table_data.append({
            "name": section["section"],
            "is_header": True
        })
        
        # Process each variable in this section
        for var in section["variables"]:
            # Get the data for each group
            if var["type"] == "continuous":
                all_data = var["format"](analyzed_fis)
                nd_data = var["format"](nd_group)
                non_nd_data = var["format"](non_nd_group)
                
                # Update the name with actual sample count
                var_name = var["name"].replace("n = *", f"n = {len(all_data)}")
                
                # Calculate result for continuous variable
                result = process_continuous_var(var_name, all_data, nd_data, non_nd_data, var["test"])
            else:  # categorical
                all_count = var["format"](analyzed_fis)
                nd_count = var["format"](nd_group)
                non_nd_count = var["format"](non_nd_group)
                
                # Update the name with actual sample count
                var_name = var["name"].replace("n = *", f"n = {total_n}")
                
                # Process categorical variable
                result = process_categorical_var(var_name, all_count, nd_count, non_nd_count, 
                                               total_n, nd_n, non_nd_n, var["test"])
            
            table_data.append(result)
    
    # Create the table
    create_table_pdf(table_data, total_n, nd_n, non_nd_n, os.path.join(tables_dir, file_name))
    
    print(f"Demographics table saved to {os.path.join(tables_dir, file_name)}")
    return table_data

def process_continuous_var(name, all_data, nd_data, non_nd_data, test_type):
    """Process a continuous variable and decide appropriate test"""
    result = {"name": name, "is_header": False}
    p_value = None
    p_value_float = None
    test_symbol = ""

    # If no data for either group, return empty result
    if not nd_data or not non_nd_data:
        return {
            "name": name,
            "is_header": False,
            "all": "Insufficient data",
            "nd": "Insufficient data",
            "non_nd": "Insufficient data",
            "p_value": "-",
            "p_value_float": None, # Add float value
            "test_symbol": ""
        }

    # Decide which test to use if auto
    if test_type == "auto":
        if len(nd_data) >= 3 and len(non_nd_data) >= 3:
            _, p_nd = stats.shapiro(nd_data)
            _, p_non_nd = stats.shapiro(non_nd_data)
            if p_nd > 0.05 and p_non_nd > 0.05:
                test_type = "ttest"
            else:
                test_type = "mwu"
        else:
            test_type = "mwu"

    # Get common stats
    mean_all = np.mean(all_data)
    median_all = np.median(all_data)
    q1_all = np.percentile(all_data, 25)
    q3_all = np.percentile(all_data, 75)
    mean_nd = np.mean(nd_data)
    median_nd = np.median(nd_data)
    q1_nd = np.percentile(nd_data, 25)
    q3_nd = np.percentile(nd_data, 75)
    mean_non_nd = np.mean(non_nd_data)
    median_non_nd = np.median(non_nd_data)
    q1_non_nd = np.percentile(non_nd_data, 25)
    q3_non_nd = np.percentile(non_nd_data, 75)

    # Run the appropriate test
    if test_type == "ttest":
        _, p_var = stats.levene(nd_data, non_nd_data)
        equal_var = p_var >= 0.05
        t_stat, p_value = stats.ttest_ind(nd_data, non_nd_data, equal_var=equal_var)
        p_value_float = p_value # Store raw float
        test_symbol = "*"

        sem_all = np.std(all_data, ddof=1) / np.sqrt(len(all_data)) if len(all_data) > 0 else 0
        sem_nd = np.std(nd_data, ddof=1) / np.sqrt(len(nd_data)) if len(nd_data) > 0 else 0
        sem_non_nd = np.std(non_nd_data, ddof=1) / np.sqrt(len(non_nd_data)) if len(non_nd_data) > 0 else 0

        result["all"] = f"{mean_all:.1f} ({sem_all:.1f})"
        result["nd"] = f"{mean_nd:.1f} ({sem_nd:.1f})"
        result["non_nd"] = f"{mean_non_nd:.1f} ({sem_non_nd:.1f})"

    else:  # MWU
        try:
             # Ensure alternative and method are correctly specified based on scipy version if needed
             u_stat, p_value = stats.mannwhitneyu(nd_data, non_nd_data, alternative='two-sided') # method="auto" is default or may raise error
             p_value_float = p_value # Store raw float
        except ValueError: # Handle cases where one sample is constant
            p_value = 1.0 # Or np.nan, depending on desired handling
            p_value_float = p_value
        test_symbol = "†"

        result["all"] = f"{median_all:.1f} [{q1_all:.1f}, {q3_all:.1f}]"
        result["nd"] = f"{median_nd:.1f} [{q1_nd:.1f}, {q3_nd:.1f}]"
        result["non_nd"] = f"{median_non_nd:.1f} [{q1_non_nd:.1f}, {q3_non_nd:.1f}]"

    # Final formatting
    result["p_value"] = f"{p_value:.3f}" if p_value is not None else "-"
    result["p_value_float"] = p_value_float
    result["test_symbol"] = test_symbol
    result["n_all"] = len(all_data)
    result["n_nd"] = len(nd_data)
    result["n_non_nd"] = len(non_nd_data)

    return result

def process_categorical_var(name, all_count, nd_count, non_nd_count, total_n, nd_n, non_nd_n, test_type):
    """Process a categorical variable and decide appropriate test"""
    p_value = None
    p_value_float = None
    test_symbol = ""

    # Calculate percentages
    all_percentage = (all_count / total_n) * 100 if total_n > 0 else 0
    nd_percentage = (nd_count / nd_n) * 100 if nd_n > 0 else 0
    non_nd_percentage = (non_nd_count / non_nd_n) * 100 if non_nd_n > 0 else 0

    # Create contingency table only if counts are valid
    if nd_n > 0 and non_nd_n > 0:
        contingency_table = np.array([
            [nd_count, nd_n - nd_count],
            [non_nd_count, non_nd_n - non_nd_count]
        ])

        # Select appropriate test
        if test_type == "auto" or test_type == "chi2":
            row_totals = contingency_table.sum(axis=1)
            col_totals = contingency_table.sum(axis=0)
            # Check for zero sums before calculating expected
            if np.all(row_totals > 0) and np.all(col_totals > 0):
                 expected = np.outer(row_totals, col_totals) / contingency_table.sum()
                 if np.any(expected < 5):
                     _, p_value = stats.fisher_exact(contingency_table)
                     test_symbol = "§"
                 else:
                     # Use correction=True for Yates' correction if desired
                     chi2, p_value, _, _ = stats.chi2_contingency(contingency_table, correction=False)
                     test_symbol = "‡"
            else: # Handle cases with zero rows/columns if they occur
                 p_value = 1.0 # Or np.nan
                 test_symbol = "" # No test applicable
        else: # Fisher's exact test explicitly specified
            _, p_value = stats.fisher_exact(contingency_table)
            test_symbol = "§"

        p_value_float = p_value # Store raw float
    else:
        # Cannot perform test if one group has zero count
        p_value = None
        p_value_float = None
        test_symbol = ""


    # Format the result
    result = {
        "name": name,
        "is_header": False,
        "all": f"{all_count}/{total_n} ({all_percentage:.1f}%)",
        "nd": f"{nd_count}/{nd_n} ({nd_percentage:.1f}%)" if nd_n > 0 else "0/0 (N/A)",
        "non_nd": f"{non_nd_count}/{non_nd_n} ({non_nd_percentage:.1f}%)" if non_nd_n > 0 else "0/0 (N/A)",
        "p_value": f"{p_value:.3f}" if p_value is not None else "-",
        "p_value_float": p_value_float, # Add float value
        "test_symbol": test_symbol,
        # Store actual counts for reference if needed
        "count_all": all_count,
        "count_nd": nd_count,
        "count_non_nd": non_nd_count
    }

    return result



def process_continuous_var(name, all_data, nd_data, non_nd_data, test_type):
    """Process a continuous variable and decide appropriate test"""
    result = {"name": name, "is_header": False}
    p_value = None
    p_value_float = None
    test_symbol = ""

    # If no data for either group, return empty result
    if not nd_data or not non_nd_data:
        return {
            "name": name,
            "is_header": False,
            "all": "Insufficient data",
            "nd": "Insufficient data",
            "non_nd": "Insufficient data",
            "p_value": "-",
            "p_value_float": None, # Add float value
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
            except ValueError: # Handle potential errors in shapiro test
                test_type = "mwu"
        else:
            test_type = "mwu" # Default to non-parametric for small or constant samples

    # Get common stats only if data exists
    mean_all, median_all, q1_all, q3_all = (np.nan,) * 4
    if all_data:
        mean_all = np.mean(all_data)
        median_all = np.median(all_data)
        q1_all = np.percentile(all_data, 25)
        q3_all = np.percentile(all_data, 75)

    mean_nd, median_nd, q1_nd, q3_nd = (np.nan,) * 4
    if nd_data:
        mean_nd = np.mean(nd_data)
        median_nd = np.median(nd_data)
        q1_nd = np.percentile(nd_data, 25)
        q3_nd = np.percentile(nd_data, 75)

    mean_non_nd, median_non_nd, q1_non_nd, q3_non_nd = (np.nan,) * 4
    if non_nd_data:
        mean_non_nd = np.mean(non_nd_data)
        median_non_nd = np.median(non_nd_data)
        q1_non_nd = np.percentile(non_nd_data, 25)
        q3_non_nd = np.percentile(non_nd_data, 75)


    # Run the appropriate test
    if test_type == "ttest":
        try:
            _, p_var = stats.levene(nd_data, non_nd_data)
            equal_var = p_var >= 0.05
            t_stat, p_value = stats.ttest_ind(nd_data, non_nd_data, equal_var=equal_var, nan_policy='omit')
            p_value_float = p_value # Store raw float
            test_symbol = "*"

            sem_all = np.std(all_data, ddof=1) / np.sqrt(len(all_data)) if len(all_data) > 0 else 0
            sem_nd = np.std(nd_data, ddof=1) / np.sqrt(len(nd_data)) if len(nd_data) > 0 else 0
            sem_non_nd = np.std(non_nd_data, ddof=1) / np.sqrt(len(non_nd_data)) if len(non_nd_data) > 0 else 0

            result["all"] = f"{mean_all:.1f} ({sem_all:.1f})" if not np.isnan(mean_all) else "N/A"
            result["nd"] = f"{mean_nd:.1f} ({sem_nd:.1f})" if not np.isnan(mean_nd) else "N/A"
            result["non_nd"] = f"{mean_non_nd:.1f} ({sem_non_nd:.1f})" if not np.isnan(mean_non_nd) else "N/A"
        except ValueError: # Handle cases like empty arrays after nan removal
             p_value = np.nan
             p_value_float = np.nan
             test_symbol = "*" # Still indicate t-test was attempted
             result["all"] = f"{mean_all:.1f} ({sem_all:.1f})" if not np.isnan(mean_all) else "N/A"
             result["nd"] = f"{mean_nd:.1f} ({sem_nd:.1f})" if not np.isnan(mean_nd) else "N/A"
             result["non_nd"] = f"{mean_non_nd:.1f} ({sem_non_nd:.1f})" if not np.isnan(mean_non_nd) else "N/A"


    else:  # MWU
        # Check if both arrays are non-empty after potential nan removal if applicable
        nd_data_valid = [x for x in nd_data if x is not None and not np.isnan(x)]
        non_nd_data_valid = [x for x in non_nd_data if x is not None and not np.isnan(x)]

        if nd_data_valid and non_nd_data_valid:
            try:
                 u_stat, p_value = stats.mannwhitneyu(nd_data_valid, non_nd_data_valid, alternative='two-sided')
                 p_value_float = p_value # Store raw float
            except ValueError: # Handle cases where one sample is constant or other MWU errors
                p_value = np.nan
                p_value_float = np.nan
        else: # Not enough data for test
            p_value = np.nan
            p_value_float = np.nan

        test_symbol = "†"

        result["all"] = f"{median_all:.1f} [{q1_all:.1f}, {q3_all:.1f}]" if not np.isnan(median_all) else "N/A"
        result["nd"] = f"{median_nd:.1f} [{q1_nd:.1f}, {q3_nd:.1f}]" if not np.isnan(median_nd) else "N/A"
        result["non_nd"] = f"{median_non_nd:.1f} [{q1_non_nd:.1f}, {q3_non_nd:.1f}]" if not np.isnan(median_non_nd) else "N/A"

    # Final formatting
    result["p_value"] = f"{p_value:.3f}" if p_value is not None and not np.isnan(p_value) else "-"
    result["p_value_float"] = p_value_float
    result["test_symbol"] = test_symbol
    result["n_all"] = len(all_data)
    result["n_nd"] = len(nd_data)
    result["n_non_nd"] = len(non_nd_data)

    return result

def process_categorical_var(name, all_count, nd_count, non_nd_count, total_n, nd_n, non_nd_n, test_type):
    """Process a categorical variable and decide appropriate test"""
    p_value = None
    p_value_float = None
    test_symbol = ""

    # Calculate percentages
    all_percentage = (all_count / total_n) * 100 if total_n > 0 else 0
    nd_percentage = (nd_count / nd_n) * 100 if nd_n > 0 else 0
    non_nd_percentage = (non_nd_count / non_nd_n) * 100 if non_nd_n > 0 else 0

    # Create contingency table only if counts are valid
    # Ensure non-negative counts for subtraction
    nd_n_pos = max(0, nd_n)
    non_nd_n_pos = max(0, non_nd_n)
    nd_count_pos = max(0, nd_count)
    non_nd_count_pos = max(0, non_nd_count)

    if nd_n_pos > 0 and non_nd_n_pos > 0:
        contingency_table = np.array([
            [nd_count_pos, nd_n_pos - nd_count_pos],
            [non_nd_count_pos, non_nd_n_pos - non_nd_count_pos]
        ])

        # Select appropriate test only if table is valid (e.g., non-zero entries if required by test)
        if contingency_table.sum() > 0 and np.all(contingency_table >= 0):
            try:
                if test_type == "auto" or test_type == "chi2":
                    row_totals = contingency_table.sum(axis=1)
                    col_totals = contingency_table.sum(axis=0)
                    if np.all(row_totals > 0) and np.all(col_totals > 0):
                         expected = np.outer(row_totals, col_totals) / contingency_table.sum()
                         if np.any(expected < 5):
                             _, p_value = stats.fisher_exact(contingency_table)
                             test_symbol = "§"
                         else:
                             chi2, p_value, _, _ = stats.chi2_contingency(contingency_table, correction=False)
                             test_symbol = "‡"
                    else: # Handle cases with zero rows/columns
                         p_value = np.nan
                         test_symbol = "" # No test applicable
                else: # Fisher's exact test explicitly specified
                    _, p_value = stats.fisher_exact(contingency_table)
                    test_symbol = "§"

                p_value_float = p_value # Store raw float

            except ValueError: # Catch potential errors in stat tests
                p_value = np.nan
                p_value_float = np.nan
                test_symbol = "" # Indicate test failed or not applicable

        else: # Cannot perform test if table sum is zero or invalid
            p_value = np.nan
            p_value_float = np.nan
            test_symbol = ""
    else:
        # Cannot perform test if one group has zero count
        p_value = np.nan
        p_value_float = np.nan
        test_symbol = ""


    # Format the result
    result = {
        "name": name,
        "is_header": False,
        "all": f"{all_count}/{total_n} ({all_percentage:.1f}%)",
        "nd": f"{nd_count}/{nd_n} ({nd_percentage:.1f}%)" if nd_n > 0 else "0/0 (N/A)",
        "non_nd": f"{non_nd_count}/{non_nd_n} ({non_nd_percentage:.1f}%)" if non_nd_n > 0 else "0/0 (N/A)",
        "p_value": f"{p_value:.3f}" if p_value is not None and not np.isnan(p_value) else "-",
        "p_value_float": p_value_float, # Add float value
        "test_symbol": test_symbol,
        # Store actual counts for reference if needed
        "count_all": all_count,
        "count_nd": nd_count,
        "count_non_nd": non_nd_count
    }

    return result


def create_table_pdf(table_data, total_n, nd_n, non_nd_n, output_path):
    """Creates and saves a formatted demographics table PDF."""
    columns = ["Variable", f"Entire sample (n = {total_n})", f"ND (n = {nd_n})", f"No ND (n = {non_nd_n})", "p-value"]

    # --- Data Preparation ---
    df_display_rows = []
    p_value_metadata = []
    section_header_indices_table = []
    current_table_row_index = 1
    first_data_header_skipped = False
    original_data_index = 0

    while original_data_index < len(table_data):
        row_data = table_data[original_data_index]
        if row_data.get("is_header", False):
            if not first_data_header_skipped:
                first_data_header_skipped = True
                original_data_index += 1
                continue
            df_display_rows.append([row_data["name"], "", "", "", ""])
            section_header_indices_table.append(current_table_row_index)
            p_value_metadata.append(None)
        else:
            p_value_str = row_data["p_value"]
            p_value_float = row_data.get("p_value_float", None)
            test_symbol = row_data.get('test_symbol', '')
            df_display_rows.append([
                row_data["name"], row_data["all"], row_data["nd"],
                row_data["non_nd"], p_value_str
            ])
            p_value_metadata.append({'float': p_value_float, 'symbol': test_symbol})
        current_table_row_index += 1
        original_data_index += 1

    # --- Table Creation ---
    base_fontsize = 8.5
    footnote_fontsize = 8
    row_height = 0.28
    header_height = 0.5
    title_space = 0.3
    footer_space = 1.2
    fig_height = len(df_display_rows) * row_height + header_height + title_space + footer_space

    fig, ax = plt.subplots(figsize=(12, fig_height))
    ax.axis('off')

    # Define subplot margins - needed for title positioning
    left_margin = 0.05
    right_margin = 0.95
    bottom_margin = 0.15 # Estimate based on footnote space
    top_margin = 0.95 # Estimate space below title

    plt.subplots_adjust(left=left_margin, right=right_margin, bottom=bottom_margin, top=top_margin)


    table = ax.table(
        cellText=df_display_rows,
        colLabels=columns,
        loc='center', # Table centered within the axes
        cellLoc='left',
        colWidths=[0.35, 0.18, 0.18, 0.18, 0.11]
    )

    table.auto_set_font_size(False)
    table.set_fontsize(base_fontsize)
    table.scale(1, 1.5)

    # --- Cell Formatting ---
    for i in range(len(df_display_rows) + 1):
        for j in range(len(columns)):
            try:
                cell = table[(i, j)]
            except KeyError:
                print(f"Warning: Could not get cell ({i}, {j})")
                continue

            cell.visible_edges = 'BTRL'
            cell.PAD = 0.02
            cell._loc = 'left'
            text_obj = cell.get_text()
            text_obj.set_horizontalalignment('left')
            text_obj.set_verticalalignment('center')

            if i == 0: # Header Row
                cell.set_text_props(weight='bold', style='italic')
            else: # Data/Section Header Rows
                row_index_in_lists = i - 1
                if i in section_header_indices_table: # Section Header
                    cell.set_text_props(weight='bold', style='italic')
                    if j != 0 and j != 4:
                         cell.get_text().set_text("")
                         cell.set_facecolor('white')
                         cell.visible_edges = 'TB'
                    elif j == 4:
                        cell.get_text().set_text("")
                        cell.set_facecolor('white')
                        cell.visible_edges = 'TBR'
                    else:
                         cell.visible_edges = 'BLT'
                else: # Data Row
                    if j == 4: # p-value column
                        meta = p_value_metadata[row_index_in_lists]
                        if meta:
                            p_val_float = meta['float']
                            test_symbol = meta['symbol']
                            p_val_str = cell.get_text().get_text()
                            is_significant = False
                            # Check significance, handling potential non-numeric p_val_str
                            if p_val_float is not None and not np.isnan(p_val_float):
                                if p_val_float < 0.001:
                                     p_val_str = "<0.001" # Use exact string for display
                                     is_significant = True
                                elif p_val_float < 0.05:
                                     is_significant = True

                            p_val_fmt = p_val_str
                            if is_significant:
                                p_val_fmt = r"\mathbf{" + p_val_fmt + "}"
                            if test_symbol:
                                formatted_symbol = test_symbol.replace('†', r'\dagger').replace('‡', r'\ddagger').replace('§', r'\S')
                                p_val_fmt += f"^{{{formatted_symbol}}}"
                            cell.get_text().set_text(f"${p_val_fmt}$")

    # --- Title ---
    title_text = r"$\bf{Table\ 1:}$ Sample demographic and clinical characteristics"
    plt.figtext(left_margin, top_margin - .07,
                title_text,
                ha='left', va='bottom',
                fontsize=base_fontsize)

    # --- Footnotes ---
    footnote_y_start = 0.07 # Original y position
    footnote_spacing = 0.025
    nb_text = "NB: percentages reported reflect the percentage of measured values falling into a given category (i.e., missing or uncoded values are not included in the denominator)."
    test_symbols_text = (
        r"$^*$Student$'$s t-test used" + "\n" +
        r"$^\dagger$Mann-Whitney U-test used" + "\n" +
        r"$^\ddagger$Chi-square test used" + "\n" +
        r"$^\S$Fisher$'$s exact test used"
    )
    # Place footnotes using original relative positions, aligned left
    # Use left_margin for x coordinate to align with title and table area
    fig.text(left_margin, footnote_y_start, nb_text, # Remove textwrap.fill
             ha='left', va='top', fontsize=footnote_fontsize)
    fig.text(left_margin, footnote_y_start - footnote_spacing, test_symbols_text,
             ha='left', va='top', fontsize=footnote_fontsize)

    # --- Layout and Save ---
    # Adjust bottom margin in subplots_adjust if footnotes are clipped
    lowest_footnote_pos = footnote_y_start - footnote_spacing - 0.04 # Estimate needed space
    plt.subplots_adjust(bottom=max(0.01, lowest_footnote_pos)) # Ensure bottom > 0

    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close(fig)

    print(f"Demographics table saved to {output_path}")



def generate_demographics_tables(stat_gui_instance, case = None):
    """
    Function to be called from stat_gui.py to generate demographics tables
    
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
    create_demographics_table(stat_gui_instance.fis_dict, base_save_dir, case = case)
    
    print("Demographics tables generated successfully.")