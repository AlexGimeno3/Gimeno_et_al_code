import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import textwrap
import re
from PyPDF2 import PdfMerger
import sys
current_dir = os.path.dirname(os.path.abspath(__file__))
utils_dir = os.path.join(current_dir, '..', '..', 'utils')
sys.path.append(utils_dir)
import tempfile
import subprocess

def get_base_save_dir(start_dir):
    """Run file dialog in completely separate Python process."""
    script_content = f'''
import tkinter as tk
from tkinter import filedialog
import os
import sys

try:
    root = tk.Tk()
    root.withdraw()
    root.lift()
    root.attributes('-topmost', True)
    
    start_dir = r"{start_dir}"
    selected = filedialog.askdirectory(
        title="Select a subdirectory",
        initialdir=start_dir if os.path.exists(start_dir) else None
    )
    
    if selected:
        print(selected)
    else:
        print("CANCELLED")
        
    root.destroy()
    
except Exception as e:
    print(f"ERROR: {{e}}", file=sys.stderr)
    sys.exit(1)
'''
    
    try:
        # Write script to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
            f.write(script_content)
            script_path = f.name
        
        # Run in separate process
        result = subprocess.run(
            [sys.executable, script_path], 
            capture_output=True, 
            text=True, 
            timeout=6000000
        )
        
        # Clean up temp file
        os.unlink(script_path)
        
        if result.returncode == 0:
            output = result.stdout.strip()
            return None if output == "CANCELLED" else output
        else:
            print(f"Subprocess error: {result.stderr}")
            return None
            
    except subprocess.TimeoutExpired:
        print("File dialog timed out")
        return None
    except Exception as e:
        print(f"Subprocess method failed: {e}")
        return None

#Sort categories function (to return variables for delta, theta, all, etc)
def sort_by_categories(df_list):
    # Extract headers and rows
    headers = df_list[0]
    rows = df_list[1:]
    # Find the column index for variable name (assuming it's called 'Variable' or similar)
    var_idx = next((i for i, col in enumerate(headers) if 'variable' in col.lower()), 0)
    p_idx = next((i for i, col in enumerate(headers) if 'test results' in col.lower()), -1)
    # Define categories
    categories = {
        'delta': [],
        'theta': [],
        'alpha': [],
        'beta': [],
        'other': []
    }
    # Categorize rows
    for row in rows:
        var_name = str(row[var_idx]).lower()
        if 'delta' in var_name:
            categories['delta'].append(row)
        elif 'theta' in var_name:
            categories['theta'].append(row)
        elif 'alpha' in var_name:
            categories['alpha'].append(row)
        elif 'beta' in var_name:
            categories['beta'].append(row)
        else:
            categories['other'].append(row)
    # Sort each category alphabetically by variable name
    # for category in categories:
    #     categories[category].sort(key=lambda x: str(x[var_idx]))
    # Combine all sorted categories
    sorted_rows = []
    for category in ['delta', 'theta', 'alpha', 'beta', 'other']:
        sorted_rows.extend(categories[category])
    return headers, sorted_rows, p_idx

# Convert to pandas DataFrame
def list_to_df(headers, rows):
    return pd.DataFrame(rows, columns=headers)

# Generate PDF table with formatting
def generate_pdf_table(df, p_idx, title, filename, custom_group = None):
    greek_dict = {
        "alpha": "α",
        "beta": "β",
        "theta": "θ",
        "tau": "τ",
        "delta": "δ"
    }
    
    # Helper function to apply replacements
    def apply_replacements(text):
        if not isinstance(text, str):
            text = str(text)
        
        # Replace reeg with rEEG (case-insensitive)
        text = re.sub(r'\breeg\b', 'rEEG', text, flags=re.IGNORECASE)
        
        # Replace Greek letters (case-insensitive)
        for greek_key, greek_val in greek_dict.items():
            text = re.sub(greek_key, greek_val, text, flags=re.IGNORECASE)
        
        return text
    
    # Apply replacements to title
    title = apply_replacements(title)
    
    # Create folder if it doesn't exist
    folder = os.path.dirname(filename)
    if folder and not os.path.exists(folder):
        os.makedirs(folder)
    
    # Group by the variable name - extract the variable part before the colon
    var_pattern = re.compile(r'^(.*?)(?:,|$)')
    if custom_group is not None:
        if custom_group == "all":
            var_pattern = re.compile(r'^(.*)$')
    
    # Create a dictionary to group rows by variable name
    var_groups = {}
    for i, row in df.iterrows():
        var_name = row['Variable']
        # Fix random colons by replacing ", :" with ","
        var_name = re.sub(r',\s*:', ',', var_name)
        
        # Apply replacements to variable names
        var_name = apply_replacements(var_name)
        
        # Store the cleaned variable name back to the dataframe
        df.at[i, 'Variable'] = var_name
        
        match = var_pattern.match(var_name)
        if match:
            base_var = match.group(1).strip()  # Get the variable name part
        else:
            base_var = var_name  # Use the whole string if no match
            
        if base_var not in var_groups:
            var_groups[base_var] = []
        var_groups[base_var].append(row)
    
    # Create a PDF with one page per variable
    with PdfPages(filename) as pdf:
        for var_name, rows in var_groups.items():
            # Create a dataframe for this variable only
            var_df = pd.DataFrame(rows)
            
            # Format numerical values in the dataframe to prevent overflow
            for col in var_df.columns:
                if col != 'Variable' and col != 'Equal Var':
                    var_df[col] = var_df[col].apply(lambda x: 
                        re.sub(r'(\d+\.\d+)', lambda match: f"{float(match.group(0)):.3f}", str(x))  # Round to 3 decimals
                        if isinstance(x, str) else x)
            
            # Set up the figure for this variable
            fig_width = min(16, 2 + len(var_df.columns) * 2)  # Adjusted width calculation
            fig_height = min(12, 1 + len(var_df) * 0.8)
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            
            ax.axis('off')
            ax.axis('tight')
            
            # Wrap text for better readability
            wrapped_text = []
            for val in var_df['Variable']:
                # Fix random colons again to be safe
                val = re.sub(r',\s*:', ',', str(val))
                wrapped = "\n".join(textwrap.wrap(str(val), width=25))
                wrapped_text.append(wrapped)
            var_df = var_df.copy()
            var_df['Variable'] = wrapped_text
            
            # Process all cells to ensure they fit - more aggressive for Theil-Sen
            for col in var_df.columns:
                if col not in ['Variable', 'Equal Var']:
                    var_df[col] = var_df[col].apply(lambda x: 
                        # Round Theil-Sen values to 3 decimal places 
                        re.sub(r'(\d+\.\d+)', lambda match: f"{float(match.group(0)):.3f}", str(x))
                        if 'theil_sen' in str(x).lower() else 
                        re.sub(r'(\d+\.\d+)', lambda match: f"{float(match.group(0)):.3f}", str(x)))  # 3 decimal places for all values

            
            # Adjust cell colors
            cell_colors = np.full((len(var_df), len(var_df.columns)), 'white', dtype=object)
            if p_idx >= 0:
                for i, row in enumerate(var_df.itertuples(index=False)):
                    p_val_str = str(row[p_idx])
                    match = re.search(r'p = (\d+\.\d+)', p_val_str)
                    if match:
                        p_val = float(match.group(1))
                    else:
                        try:
                            p_val = float(p_val_str)
                        except:
                            p_val = 1.0
                    if p_val < 0.05:
                        cell_colors[i] = 'lightgreen'
                    elif p_val < 0.10:
                        cell_colors[i] = 'lightyellow'
            
            table = ax.table(
                cellText=var_df.values,
                colLabels=var_df.columns,
                cellColours=cell_colors,
                loc='center',
                cellLoc='center'
            )
            
            table.auto_set_font_size(False)
            table.set_fontsize(7)  # Smaller font
            table.scale(1.2, 2.2)  # Adjusted scaling for better balance
            
            # Adjust column widths based on content; NB: values here are relative
            for (i, j), cell in table.get_celld().items():
                if j == 0:  # Variable column
                    cell.set_width(0.35)  #  variable column
                elif j == 4: #results column
                    cell.set_width(0.40)
                    cell.set_fontsize(7)
                else:
                    cell.set_width(0.27)  #   data columns
            
            # Add variable name as page title (without "All Statistical Results")
            # Apply replacements to the page title as well
            page_title = apply_replacements(var_name)
            plt.title(f"{page_title}", fontsize=14, pad=20)  # Removed title parameter
            plt.tight_layout()
            
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            
    print(f"PDF generated with {len(var_groups)} variable pages: {filename}")

# Split into categories and save into separate PDFs
def save_by_category(df_list, base_folder, base_name, custom_group = None):
    headers, sorted_rows, p_idx = sort_by_categories(df_list)
    
    # Split by categories
    categories = ['delta', 'theta', 'alpha', 'beta', 'other']
    var_idx = next((i for i, col in enumerate(headers) if 'variable' in col.lower()), 0)
    
    for cat in categories:
        if cat == 'other':
            cat_rows = [row for row in sorted_rows if not any(band in str(row[var_idx]).lower() for band in ['delta', 'theta', 'alpha', 'beta'])]
        else:
            cat_rows = [row for row in sorted_rows if cat in str(row[var_idx]).lower()]
        if cat_rows:
            cat_df = list_to_df(headers, cat_rows)
            generate_pdf_table(cat_df, p_idx, f"{cat.capitalize()} {base_name}", os.path.join(base_folder, f"{cat}_results_table.pdf"), custom_group=custom_group)

def combine_pdfs(folder, output_filename):
    merger = PdfMerger()
    pdf_files = [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith('.pdf')]
    pdf_files.sort()  # optional: alphabetical order (delta, then theta, etc.)
    for pdf in pdf_files:
        merger.append(pdf)
    merger.write(os.path.join(folder, output_filename))
    merger.close()
    print(f"Combined PDF created: {os.path.join(folder, output_filename)}")

def integrate_PDF(base_save_dir, significant_df, trending_df):
    # Ensure output directory exists
    output_dir = os.path.join(base_save_dir,'significant_results')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Combine significant and trending dataframes
    sig_trend_headers = []
    sig_trend_rows = []
    
    if significant_df:
        sig_headers = significant_df[0]
        sig_rows = significant_df[1:]
        if not sig_trend_headers:
            sig_trend_headers = sig_headers
        sig_trend_rows.extend(sig_rows)
    
    if trending_df:
        trend_headers = trending_df[0]
        trend_rows = trending_df[1:]
        if not sig_trend_headers:
            sig_trend_headers = trend_headers
        sig_trend_rows.extend(trend_rows)
    
    # Group by variable name
    var_idx = next((i for i, col in enumerate(sig_trend_headers) if 'variable' in col.lower()), 0)
    p_idx = next((i for i, col in enumerate(sig_trend_headers) if 'test results' in col.lower()), -1)
    
    # Sort and organize by categories
    categories = ['delta', 'theta', 'alpha', 'beta', 'other']
    categorized_rows = {cat: [] for cat in categories}
    
    # Categorize rows
    for row in sig_trend_rows:
        var_name = str(row[var_idx]).lower()
        if 'delta' in var_name:
            categorized_rows['delta'].append(row)
        elif 'theta' in var_name:
            categorized_rows['theta'].append(row)
        elif 'alpha' in var_name:
            categorized_rows['alpha'].append(row)
        elif 'beta' in var_name:
            categorized_rows['beta'].append(row)
        else:
            categorized_rows['other'].append(row)
    
    # Generate PDF with combined data
    combined_rows = []
    for cat in categories:
        combined_rows.extend(categorized_rows[cat])
    
    combined_df = pd.DataFrame(combined_rows, columns=sig_trend_headers)
    
    # Group data by variable base name for page generation
    import re
    var_pattern = re.compile(r'^(.*?)(?:,|$)')
    var_groups = {}
    
    for i, row in combined_df.iterrows():
        var_name = row['Variable']
        # Fix random colons by replacing ", :" with ","
        var_name = re.sub(r',\s*:', ',', var_name)
        # Store the cleaned variable name back to the dataframe
        combined_df.at[i, 'Variable'] = var_name
        
        match = var_pattern.match(var_name)
        if match:
            base_var = match.group(1).strip()  # Get the variable name part
        else:
            base_var = var_name  # Use the whole string if no match
            
        if base_var not in var_groups:
            var_groups[base_var] = []
        var_groups[base_var].append(row)
    
    # Create a PDF with one page per variable
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import numpy as np
    import textwrap
    
    output_path = os.path.join(output_dir, 'significant_and_trending_combined.pdf')
    
    with PdfPages(output_path) as pdf:
        for var_name, rows in var_groups.items():
            # Create a dataframe for this variable only
            var_df = pd.DataFrame(rows)
            
            # Format numerical values in the dataframe to prevent overflow
            for col in var_df.columns:
                if col != 'Variable' and col != 'Equal Var':
                    var_df[col] = var_df[col].apply(lambda x: 
                        re.sub(r'(\d+\.\d+)', lambda match: f"{float(match.group(0)):.3f}", str(x))
                        if isinstance(x, str) else x)
            
            # Set up the figure for this variable
            fig_width = min(18, 2 + len(var_df.columns) * 2.5)
            fig_height = min(12, 1 + len(var_df) * 0.8)
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            
            ax.axis('off')
            ax.axis('tight')
            
            # Wrap text for better readability
            wrapped_text = []
            for val in var_df['Variable']:
                val = re.sub(r',\s*:', ',', str(val))
                wrapped = "\n".join(textwrap.wrap(str(val), width=30))
                wrapped_text.append(wrapped)
            var_df = var_df.copy()
            var_df['Variable'] = wrapped_text
            
            # Process all cells to ensure they fit
            for col in var_df.columns:
                if col not in ['Variable', 'Equal Var']:
                    var_df[col] = var_df[col].apply(lambda x: 
                        re.sub(r'(\d+\.\d+)', lambda match: f"{float(match.group(0)):.3f}", str(x)))
            
            # Adjust cell colors based on p-value
            cell_colors = np.full((len(var_df), len(var_df.columns)), 'white', dtype=object)
            if p_idx >= 0:
                for i, row in enumerate(var_df.itertuples(index=False)):
                    p_val_str = str(row[p_idx])
                    match = re.search(r'p = (\d+\.\d+)', p_val_str)
                    if match:
                        p_val = float(match.group(1))
                    else:
                        try:
                            p_val = float(p_val_str)
                        except:
                            p_val = 1.0
                    if p_val < 0.05:
                        cell_colors[i] = 'lightgreen'  # Significant
                    elif p_val < 0.10:
                        cell_colors[i] = 'lightyellow'  # Trending
            
            table = ax.table(
                cellText=var_df.values,
                colLabels=var_df.columns,
                cellColours=cell_colors,
                loc='center',
                cellLoc='center'
            )
            
            table.auto_set_font_size(False)
            table.set_fontsize(9)
            table.scale(1.3, 2.5)
            
            # Adjust column widths based on content
            for (i, j), cell in table.get_celld().items():
                if j == 0:  # Variable column
                    cell.set_width(0.35)
                else:
                    cell.set_width(0.15)
            
            # Add variable name as page title
            plt.title(f"{var_name}", fontsize=14, pad=20)
            plt.tight_layout()
            
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
    
    print(f"Integrated significant and trending results PDF created: {output_path}")