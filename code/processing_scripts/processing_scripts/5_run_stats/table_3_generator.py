import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


def create_seizure_table(fis_dict, base_save_dir, file_name="table_3_seizure_data.pdf"):
    """
    Create and save a table showing seizure data by neurodevelopmental outcome.
    
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
    analyzed_fis_clean = [fis_obj for fis_obj in fis_dict.values() if fis_obj.analyzed_whole_sample == 1 and not fis_obj.seizures_intraop is None]
    analyzed_fis = [fis_obj for fis_obj in fis_dict.values() if fis_obj.analyzed_whole_sample == 1]
    
    
    if not analyzed_fis:
        print("No analyzed FIS objects found. Cannot create seizure table.")
        return
    
    # Split into ND and non-ND groups
    nd_group = [fis_obj for fis_obj in analyzed_fis if fis_obj.ND == 1]
    non_nd_group = [fis_obj for fis_obj in analyzed_fis if fis_obj.ND == 0]
    
    # Count seizures in each group
    # Combining intraop and postop seizures for this table
    nd_seizure_count = sum(1 for fis_obj in nd_group if fis_obj.seizures_intraop == 1)
    nd_no_seizure_count = sum(1 for fis_obj in nd_group if fis_obj.seizures_intraop == 0)
    non_nd_seizure_count = sum(1 for fis_obj in non_nd_group if fis_obj.seizures_intraop == 1)
    non_nd_no_seizure_count = sum(1 for fis_obj in non_nd_group if fis_obj.seizures_intraop == 0)
    
    # Create contingency table for Fisher's exact test
    contingency_table = np.array([
        [nd_seizure_count, nd_no_seizure_count],
        [non_nd_seizure_count, non_nd_no_seizure_count]
    ])
    
    # Perform Fisher's exact test
    _, p_value = stats.fisher_exact(contingency_table)
    
    # Generate the PDF table
    create_seizure_table_pdf(
        nd_n=len(nd_group),
        non_nd_n=len(non_nd_group),
        nd_seizure_count=nd_seizure_count,
        nd_no_seizure_count=nd_no_seizure_count,
        non_nd_seizure_count=non_nd_seizure_count,
        non_nd_no_seizure_count=non_nd_no_seizure_count,
        p_value=p_value,
        output_path=os.path.join(tables_dir, file_name)
    )
    
    print(f"Seizure table saved to {os.path.join(tables_dir, file_name)}")


def create_seizure_table_pdf(nd_n, non_nd_n, nd_seizure_count, nd_no_seizure_count, 
                             non_nd_seizure_count, non_nd_no_seizure_count, p_value, output_path):
    """
    Creates and saves a formatted seizure table PDF.
    
    Parameters:
    -----------
    nd_n : int
        Total number of subjects in ND group
    non_nd_n : int
        Total number of subjects in non-ND group
    nd_seizure_count : int
        Number of subjects with seizures in ND group
    nd_no_seizure_count : int
        Number of subjects without seizures in ND group
    non_nd_seizure_count : int
        Number of subjects with seizures in non-ND group
    non_nd_no_seizure_count : int
        Number of subjects without seizures in non-ND group
    p_value : float
        P-value from Fisher's exact test
    output_path : str
        Path to save the PDF file
    """
    # Calculate percentages
    nd_seizure_pct = (nd_seizure_count / nd_n) * 100 if nd_n > 0 else 0
    nd_no_seizure_pct = (nd_no_seizure_count / nd_n) * 100 if nd_n > 0 else 0
    non_nd_seizure_pct = (non_nd_seizure_count / non_nd_n) * 100 if non_nd_n > 0 else 0
    non_nd_no_seizure_pct = (non_nd_no_seizure_count / non_nd_n) * 100 if non_nd_n > 0 else 0
    
    # Format p-value
    if p_value < 0.001:
        p_value_str = "<0.001"
    else:
        p_value_str = f"{p_value:.3f}"
    
    # Create table data
    columns = [f"", f"ND (n = {nd_n})", f"No ND (n = {non_nd_n})"]
    
    table_data = [
        ["No seizure (confirmed)", f"{nd_no_seizure_count} ({nd_no_seizure_pct:.1f}%)", 
         f"{non_nd_no_seizure_count} ({non_nd_no_seizure_pct:.1f}%)"],
        ["Seizure (confirmed)", f"{nd_seizure_count} ({nd_seizure_pct:.1f}%)", 
         f"{non_nd_seizure_count} ({non_nd_seizure_pct:.1f}%)"],
        [r"$\chi^2$/Fisher's exact test", "", f"p = {p_value_str}"]
    ]
    
    # Create figure and table
    fig, ax = plt.subplots(figsize=(7, 3))
    ax.axis('off')
    
    # Create the table
    table = ax.table(
        cellText=table_data,
        colLabels=columns,
        loc='center',
        cellLoc='center'
    )
    
    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)
    
    # Apply styling to cells
    for i in range(len(table_data) + 1):  # +1 for header row
        for j in range(len(columns)):
            cell = table[(i, j)]
            cell.visible_edges = 'BTRL'  # Show all borders
            
            # Header row
            if i == 0:
                cell.set_text_props(weight='bold')
            
            # Test result row
            if i == 3:  # Fisher's test row
                cell.set_facecolor('0.9')  # Light gray background
    
    # Add title
    plt.figtext(0.05, 0.95, r"$\bf{Table\ 3:}$ Seizure data by neurodevelopmental outcome", 
                ha='left', fontsize=10)
    
    # Add footnote
    footnote = "NB: Fisher's exact test was also included, as one cell had an expected value count under 5."
    plt.figtext(0.05, 0.02, footnote, ha='left', fontsize=8)
    
    # Save the figure
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def generate_seizure_table(stat_gui_instance):
    """
    Function to be called from stat_gui.py to generate the seizure table
    
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
    
    # Create the seizure table
    create_seizure_table(stat_gui_instance.fis_dict, base_save_dir)
    
    print("Seizure table generated successfully.")