from stat_tests import var_dict, stat_dict
import numpy as np
import textwrap


def save_correlation_results(self, correlation_data, folder_name="correlation_results"):
    """
    Save correlation results to PDF files with information about whether values were adjusted.
    
    Parameters:
    -----------
    correlation_data : list
        List of correlation results tuples
    folder_name : str
        Name of the folder to save results to (default: 'correlation_results')
        Can include suffix to indicate if adjusted values were used
    """
    import os
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import numpy as np
    from PDF_fxns import save_by_category, integrate_PDF, sort_by_categories, list_to_df
    import seaborn as sns
    
    # Set professional plot style
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_palette("husl")
    
    # Set default font sizes for professional appearance
    plt.rcParams.update({
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.figsize': (8, 6),
        'figure.dpi': 300
    })

    base_dir = os.path.join(self.base_save_dir, folder_name)
    os.makedirs(base_dir, exist_ok=True)
    figures_dir = os.path.join(base_dir, "figures")
    os.makedirs(figures_dir, exist_ok=True)
    pdf_path = os.path.join(base_dir, "correlation_results_by_variable.pdf")
    
    # Check if we're using adjusted values (based on folder name)
    using_adjusted_values = folder_name.lower().endswith("_adjusted")
    adjustment_str = "Adjusted" if using_adjusted_values else ""

    # Define dictionaries for formatting
    greek_dict = {
        "alpha": "α",
        "beta": "β",
        "theta": "θ",
        "tau": "τ",
        "delta": "δ",
        "postop": "postoperative",
        "intraop": "intraoperative",
        "preop": "preoperative",
        "enolase": "neuron-specific enolase"
    }

    units_dict = {
        "enolase": "ng/mL",
        "S100B": "μg/L",
    }
    
    def format_label(text, type, add_units=True):
        """Returns axis label string
        - type "param" for an EEG parameter, "bio" for a biomarker
        """
        # Replace underscores with spaces
        formatted = text.replace("_", " ")
        
        # Replace instances of greek_dict keys (case-insensitive)
        for key, value in greek_dict.items():
            import re
            formatted = re.sub(rf'\b{key}\b', value, formatted, flags=re.IGNORECASE)
        
        if type == "param":
            # Replace ", " with " ("
            if "," in formatted:
                formatted = formatted.replace(", ", " (")
                formatted = formatted+")"
        
        # Replace hour formats
        formatted = formatted.replace(" 24 hr", "; 24 hours postoperative")
        formatted = formatted.replace(" 72 hr", "; 72 hours postoperative")
        
        # Capitalize first character only (unless it's a greek symbol)
        if formatted and formatted[0] not in "αβθτδ":
            formatted = formatted[0].upper() + formatted[1:]
        
        # Add units if applicable
        if add_units:
            for unit_key, unit_value in units_dict.items():
                if unit_key.lower() in text.lower():
                    formatted += f" ({unit_value})"
                    break
        
        return formatted

    def format_title(parameter, biomarker):
        biomarker_formatted = format_label(biomarker, "bio", add_units=False)
        return f"{format_label(parameter, "param", add_units=False)} vs. {biomarker_formatted[0].lower() + biomarker_formatted[1:]}"
    
    # Group rows by Variable
    grouped = {}
    all_rows = []
    sig_rows = []
    trend_rows = []

    for row, eeg_vals, bio_vals, biomarker in correlation_data:
        key = row["Variable"]
        if key not in grouped:
            grouped[key] = []
        grouped[key].append((row, eeg_vals, bio_vals, biomarker))

        row_vals = [row[k] for k in ["Variable", "Biomarker", "r", "r^2", "Test results (p-value)", "n"]]
        all_rows.append(row_vals)

        try:
            p = float(row["Test results (p-value)"])
        except:
            p = 1.0

        if p < 0.05:
            sig_rows.append(row_vals)
        elif p < 0.10:
            trend_rows.append(row_vals)

        if p < 0.05:
            eeg_vals = np.array(eeg_vals)
            bio_vals = np.array(bio_vals)
            
            # Create figure with professional styling
            fig, ax = plt.subplots(figsize=(8, 6))
            
            # Create scatter plot with professional colors
            scatter = ax.scatter(eeg_vals, bio_vals, 
                               alpha=0.6, 
                               s=60,
                               color='#2E86AB',
                               edgecolors='#1A5490',
                               linewidth=1,
                               label='Data points')

            # Add line of best fit
            try:
                m, b = np.polyfit(eeg_vals, bio_vals, 1)
                x_line = np.linspace(eeg_vals.min(), eeg_vals.max(), 100)
                y_line = m * x_line + b
                ax.plot(x_line, y_line, 
                       color='#E63946', 
                       linestyle='-', 
                       linewidth=2,
                       label='Line of best fit')
                
                # Add confidence interval (optional - for even more professional look)
                # This adds a shaded region showing uncertainty
                from scipy import stats
                predict_mean_se = lambda x: stats.t.ppf(0.975, len(eeg_vals)-2) * \
                                           np.sqrt(np.mean((bio_vals - (m*eeg_vals + b))**2) * \
                                           (1/len(eeg_vals) + (x-np.mean(eeg_vals))**2/np.sum((eeg_vals-np.mean(eeg_vals))**2)))
                margin = predict_mean_se(x_line)
                ax.fill_between(x_line, y_line - margin, y_line + margin, 
                              color='#E63946', alpha=0.1)
                
            except Exception as e:
                print(f"Linear fit failed for {row['Variable']} vs {biomarker}: {e}")

            # Format and set labels
            xlabel = format_label(row['Variable'], "param", add_units=False)
            ylabel = format_label(biomarker, "bio", add_units=True)
            
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            
            # Format title
            title_str = format_title(row['Variable'], biomarker)
            ax.set_title(title_str, 
                        fontsize=14, pad=20)
            
            # Add statistical annotation in top-left
            r_value = float(row['r'])
            p_value = float(row['Test results (p-value)'])
            
            # Format p-value display
            if p_value < 0.001:
                p_text = "p < 0.001"
            else:
                p_text = f"p = {p_value:.3f}"
            
            stats_text = f"r = {r_value:.3f}\n{p_text}"
            ax.text(0.02, 0.98, stats_text, 
                   transform=ax.transAxes,
                   verticalalignment='top',
                   bbox=dict(boxstyle='round,pad=0.5', 
                            facecolor='white', 
                            edgecolor='gray',
                            alpha=0.8),
                   fontsize=11,
                   fontweight='bold')
            
            # Improve grid appearance
            ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
            ax.set_axisbelow(True)
            
            # Add minor ticks for professional look
            ax.minorticks_on()
            ax.grid(which='minor', alpha=0.1, linestyle=':', linewidth=0.5)
            
            # Adjust legend
            legend = ax.legend(loc='lower right', 
                             frameon=True, 
                             fancybox=True, 
                             shadow=True,
                             borderpad=1,
                             framealpha=0.9)
            legend.get_frame().set_facecolor('white')
            legend.get_frame().set_edgecolor('gray')
            
            # Set spine properties for cleaner look
            for spine in ax.spines.values():
                spine.set_linewidth(1.5)
                spine.set_color('#666666')
            
            # Tight layout for better spacing
            plt.tight_layout()
            
            # Sanitize filename
            var_sanitized = "".join(c if c.isalnum() or c in ['_', '-'] else '_' for c in row['Variable'])
            biomarker_sanitized = "".join(c if c.isalnum() or c in ['_', '-'] else '_' for c in biomarker)
            filename = f"{var_sanitized}_{biomarker_sanitized}_{adjustment_str}.png"

            plt.savefig(os.path.join(figures_dir, filename), 
                       bbox_inches='tight', 
                       dpi=300,
                       facecolor='white',
                       edgecolor='none')
            plt.close()
    
    # Reset matplotlib parameters to defaults
    plt.rcParams.update(plt.rcParamsDefault)
    
    # === Part 1: Original output (adjusted for dynamic row height and wrapping) ===
    headers = ["Variable", "Biomarker", "r", "r^2", "Test results (p-value)", "n"]

    with PdfPages(pdf_path) as pdf:
        # Add a title page with information about adjustment
        plt.figure(figsize=(12, 6))
        plt.axis('off')
        title = f"Correlation Results by Variable ({adjustment_str} values)"
        plt.text(0.5, 0.5, title, fontsize=16, ha='center', va='center')
        if using_adjusted_values:
            plt.text(0.5, 0.4, "Values were adjusted using ANCOVA before correlation analysis", 
                    fontsize=12, ha='center', va='center', style='italic')
        pdf.savefig()
        plt.close()
        
        for variable, entries in grouped.items():
            data = []
            colors = []

            for idx, (row, _, _, _) in enumerate(entries):
                row_data = []
                row_colors = []
                for h in headers:
                    val = str(row[h])
                    if h in ["Variable", "Biomarker"]:
                        wrapped_val = "\n".join(textwrap.wrap(val, width=30))
                    else:
                        wrapped_val = val
                    row_data.append(wrapped_val)

                try:
                    p = float(row["Test results (p-value)"])
                except:
                    p = 1.0

                if idx == 0:
                    row_colors = ['white'] * len(headers)
                elif p < 0.05:
                    row_colors = ['lightgreen'] * len(headers)  # light green for significant
                elif p < 0.10:
                    row_colors = ['lightyellow'] * len(headers)  # light yellow for trending
                else:
                    row_colors = ['white'] * len(headers)

                data.append(row_data)
                colors.append(row_colors)

            # Estimate figure height based on number of rows
            row_height = 0.45  # Slightly taller to fit wrapped lines
            fig_height = 0.8 + row_height * len(data)
            fig, ax = plt.subplots(figsize=(13, fig_height))  # wider figure
            ax.axis('off')

            # Create the table
            table = ax.table(
                cellText=[headers] + data,
                cellColours=[['lightgrey'] * len(headers)] + colors,
                loc='center',
                cellLoc='center'
            )
            table.auto_set_font_size(False)
            table.set_fontsize(8)  # smaller font for compact layout
            table.scale(1.3, 2.0)  # larger vertical scaling

            # Bold header row
            for j in range(len(headers)):
                cell = table[0, j]
                cell.set_text_props(weight='bold')

            # Title for each variable with adjustment info
            plt.title(f"{variable} ({adjustment_str} values)", fontsize=14, pad=16)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)


    print(f"Correlation results PDF created: {pdf_path}")

    # === Part 2: Additional categorized outputs ===
    headers = ["Variable", "Biomarker", "r", "r^2", "Test results (p-value)", "n"]
    all_df = [headers] + all_rows
    sig_df = [headers] + sig_rows
    trend_df = [headers] + trend_rows

    # Include adjustment info in folder names
    save_by_category(all_df, os.path.join(base_dir, f"all_results_{adjustment_str}"), 
                    f"Correlation Results ({adjustment_str} values)", custom_group="all")
    save_by_category(sig_df, os.path.join(base_dir, f"significant_results_{adjustment_str}"), 
                    f"Significant Correlations ({adjustment_str} values)", custom_group="all")
    save_by_category(trend_df, os.path.join(base_dir, f"trending_results_{adjustment_str}"), 
                    f"Trending Correlations ({adjustment_str} values)", custom_group="all")

    integrate_PDF(base_dir, sig_df, trend_df)
    
    # Create an info file about the adjustment setting
    info_path = os.path.join(base_dir, "correlation_analysis_info.txt")
    with open(info_path, "w") as f:
        f.write(f"Correlation Analysis Information\n")
        f.write("==============================\n\n")
        f.write(f"Values used: {adjustment_str.upper()}\n")
        if using_adjusted_values:
            f.write("EEG values were adjusted using ANCOVA before correlation analysis with biomarkers.\n")
            f.write("This adjustment helps control for confounding variables and provides\n")
            f.write("correlations that may better reflect the true relationships between variables.\n")
        else:
            f.write("Raw EEG values were used for correlation analysis with biomarkers.\n")
            f.write("No adjustment for confounding variables was performed.\n")

def run_correlations_from_values(self, nd_vals, nd_fis, no_nd_vals, no_nd_fis, stat_name, variable_key):
    """
    Compute correlations between EEG summary values and biomarkers using nd_* and no_nd_* values.
    """
    from scipy.stats import pearsonr
    vars_dict = var_dict()
    stats_dict = stat_dict()

    all_results = []
    biomarker_names = [
        "preop_enolase", "preop_S100B", "postop_enolase", "postop_S100B",
        "enolase_24_hr", "S100B_24_hr", "enolase_72_hr", "S100B_72_hr"
    ]

    all_vals = nd_vals + no_nd_vals
    all_fis = nd_fis + no_nd_fis
    label = f"{vars_dict.get(variable_key, variable_key)}, {stats_dict.get(stat_name, stat_name)}"

    for biomarker in biomarker_names:
        eeg_vals = []
        bio_vals = []
        for fis, val in zip(all_fis, all_vals):
            fis_obj = self.fis_dict[fis]
            biomarker_val = getattr(fis_obj, biomarker, None)
            if biomarker_val is not None and not np.isnan(biomarker_val):
                eeg_vals.append(val)
                bio_vals.append(biomarker_val)

        if len(eeg_vals) >= 3:
            try:
                r, p = pearsonr(eeg_vals, bio_vals)
                r2 = r ** 2
            except Exception:
                r, r2, p = np.nan, np.nan, np.nan
        else:
            r, r2, p = np.nan, np.nan, np.nan

        result = {
            "Variable": label,
            "Biomarker": biomarker,
            "r": round(r, 3) if not np.isnan(r) else "-",
            "r^2": round(r2, 3) if not np.isnan(r2) else "-",
            "Test results (p-value)": round(p, 4) if not np.isnan(p) else "-",
            "n": len(eeg_vals)
        }

        all_results.append((result, eeg_vals, bio_vals, biomarker))

    #print(nd_fis)
    #print(no_nd_fis)
    return all_results

def correlations_without_filters(self, fis_dict, filt_dict, var, stat, time_range="surgery", cus_time_1=None, cus_time_2=None, 
                              comparison_category=["ND", "ND"], stat_name=None, use_adjusted_values=False, control_vars=None):
    """
    Compute correlations between EEG summary values and biomarkers, with option to use adjusted values.
    
    Parameters:
    -----------
    var : str
        Variable name to analyze
    stat : function
        Statistical function to apply to the time series
    time_range : str
        Time period to analyze ('surgery', 'ecc', 'clamp', 'half', 'custom')
    cus_time_1 : float
        Custom start time if time_range is 'custom'
    cus_time_2 : float
        Custom end time if time_range is 'custom'
    comparison_category : list
        Two-element list specifying the attribute to use for grouping and its display name
    stat_name : str
        Name of the statistic for display purposes
    use_adjusted_values : bool
        Whether to adjust values based on control_vars before correlation
    control_vars : dict
        Dictionary of control variables to adjust for (e.g., {'seizure': True, 'ga_eeg': True})
    """
    from copy import deepcopy
    
    nd_group = []
    no_nd_group = []
    nd_fis = []
    no_nd_fis = []
    nd_ctrl_values = []
    no_nd_ctrl_values = [] 
    
    for fis_obj in fis_dict.values():
        # Reset filters that should NOT preclude analysis
        temp_fis = fis_obj.fis
        filt_dict.set_filter(temp_fis, "ND", 1)
        filt_dict.set_filter(temp_fis, "seizure", 1)
        
        if filt_dict.get_filter(temp_fis) == 1:
            # Set time range based on selected option
            if time_range == "surgery" or time_range == "all":
                start_time = fis_obj.cx_start_time
                end_time = fis_obj.cx_end_time
            elif time_range == "ecc":
                if not (fis_obj.ECC_times["ecc_start"] is None or fis_obj.ECC_times["ecc_end"] is None):
                    start_time = fis_obj.ECC_times["ecc_start"]
                    end_time = fis_obj.ECC_times["ecc_end"]
            elif time_range == "half":
                if not (fis_obj.ECC_times["ecc_start"] is None or fis_obj.ECC_times["ecc_end"] is None):
                    start_time = fis_obj.ECC_times["ecc_start"]
                    end_time = fis_obj.ECC_times["ecc_end"]
                    end_time = start_time + (end_time - start_time) / 2  # the midpoint between start_time and end_time
            elif time_range == "clamp":
                if not (fis_obj.ECC_times["clamp_start"] is None or fis_obj.ECC_times["clamp_end"] is None):
                    start_time = fis_obj.ECC_times["clamp_start"]
                    end_time = fis_obj.ECC_times["clamp_end"]
            elif time_range == "custom":
                start_time = cus_time_1
                end_time = cus_time_2
                
            att_str = comparison_category[0]
            head_str = comparison_category[1]
            
            # Retrieve time series data 
            var_times, var_vals = fis_obj.ts_dict[var].get_vals(start_time, end_time, fis_obj.op_status)
            add_stat = stat(var_vals, var_times)
            
            # Dynamically retrieve the attribute using att_str
            att_val = getattr(fis_obj, att_str, None)
            
            # Collect control variable values for adjustment if needed
            if use_adjusted_values and control_vars and any(control_vars.values()):
                ctrl_values = {}
                for var_name, enabled in control_vars.items():
                    if enabled:
                        ctrl_values[var_name] = getattr(fis_obj, var_name, None)
                
                # Group the data based on attribute value
                if att_val == 1 and att_val is not None:
                    nd_group.append(add_stat)
                    nd_fis.append(fis_obj.fis)
                    nd_ctrl_values.append(ctrl_values)
                else:
                    no_nd_group.append(add_stat)
                    no_nd_fis.append(fis_obj.fis)
                    no_nd_ctrl_values.append(ctrl_values)
            else:
                # For all cases, add to no_nd_group for correlation analysis
                # (we'll use all data regardless of ND status for correlations)
                no_nd_group.append(add_stat)
                no_nd_fis.append(fis_obj.fis)
    
    # If using adjusted values and we have control variables
    if use_adjusted_values and control_vars and any(control_vars.values()):
        from stat_gui import adjust_for_seizure, adjust_for_continuous_variables
        
        # Apply adjustments based on selected control variables
        if control_vars.get('seizures_intraop', False):
            nd_group, no_nd_group = adjust_for_seizure(nd_group, no_nd_group, nd_ctrl_values, no_nd_ctrl_values)
        
        # Create a dict for continuous variable adjustments
        continuous_adjustments = {
            'ga_eeg': control_vars.get('ga_eeg', False),
            'pma_eeg': control_vars.get('pma_eeg', False)
        }
        
        # Only run if any continuous adjustments are selected
        if any(continuous_adjustments.values()):
            nd_group, no_nd_group = adjust_for_continuous_variables(
                nd_group, no_nd_group, nd_ctrl_values, no_nd_ctrl_values,
                control_for_ga=continuous_adjustments['ga_eeg'],
                control_for_pma=continuous_adjustments['pma_eeg']
            )
    
    # Combine all values for correlation - these are now adjusted if use_adjusted_values=True
    all_vals = nd_group + no_nd_group
    all_fis = nd_fis + no_nd_fis
    
    return run_correlations_from_values(self, nd_group, nd_fis, no_nd_group, no_nd_fis, stat_name, var)
