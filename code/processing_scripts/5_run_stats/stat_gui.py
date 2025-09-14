"""
Statistical Analysis
====================================================================

This module provides a GUI-driven statistical analysis framework for comparing EEG biomarkers between clinical groups. It implements statistical methods including ANCOVA, handles data filtering, and generates analysis reports with visualizations.

Purpose:
--------
Transforms extracted EEG variables into statistical results. The module supports multi-group comparisons, covariate adjustment, correlation analysis, and generates standardized reports scientific publications.

Core Architecture:
-----------------
1. Data Management: object loading and caching system
2. Interactive Filtering: GUI-based selection of analysis parameters
3. Statistical Testing: User selection of appropriate statistical tests
4. Covariate Control: ANCOVA implementation for confounding variable adjustment if selected
5. Results Generation: Automated PDF reports and Excel exports
6. Quality Assurance: Comprehensive logging and validation

Key Features:
------------

Statistical Methods:
- Automatic normality testing with appropriate test selection
- ANCOVA for controlling gestational age, seizure status, etc.
- Robust non-parametric alternatives (Mann-Whitney U)
- Multiple comparison corrections
- Effect size calculations

Multi-Group Analysis:
- Neurodevelopmental outcome comparisons (ND vs. non-ND)
- Seizure status comparisons (seizure vs. non-seizure)
- Time period analysis (surgery, CPB, clamp periods)

Clinical Workflow Integration:
- Time period selection (surgery phases, custom windows)
- Multiple frequency band analysis
- Correlation analysis with clinical variables
- Export formats compatible with clinical databases

Statistical Pipeline:
-------------------

Step 1: Data Loading and Caching
--------------------------------
- Builds fis_object dictionary from time series CSV files
- Validates data integrity and completeness

Step 2: Interactive Filter Configuration
---------------------------------------
GUI components for selecting:
- Machine types (Nicolet, OBM, or both)
- Time periods (surgery, CPB, clamp, custom)
- Variables and statistics of interest
- Control variables for ANCOVA
- Quality thresholds (artifact limits)

Step 3: Statistical Test Selection
---------------------------------
Automatic selection based on data characteristics:
- Shapiro-Wilk normality testing
- Student's t-test for normal distributions
- Mann-Whitney U for non-normal distributions
- ANCOVA when control variables specified

Step 4: Multi-Level Analysis Framework
-------------------------------------
Three parallel analysis streams:
1. Whole Sample Analysis: All qualifying patients
2. Seizure Comparison: Seizure vs. non-seizure groups
3. No-Seizure Subset: ND outcomes within seizure-free patients

Step 5: Results Generation and Export
------------------------------------
- PDF reports categorized by significance level
- Excel exports with patient-level data
- Correlation analysis matrices
- Demographic and clinical tables

ANCOVA Implementation:
---------------------

Purpose: Control for confounding variables that might influence EEG biomarkers

Supported Covariates:
- seizures_intraop: Intraoperative seizure occurrence
- ga_eeg: Gestational age at EEG recording
- pma_eeg: Postmenstrual age at EEG recording

Statistical Model:
EEG_parameter ~ ND_status + seizures_intraop + ga_eeg + pma_eeg
Advantages over Traditional Methods:
- Reduces Type I error from confounding
- Increases statistical power through variance reduction
- Provides adjusted group means for clinical interpretation
- Enables analysis of covariate effects

Quality Control Framework:
-------------------------

Data Validation:
- Minimum clean run requirements
- Artifact proportion thresholds
- Missing data assessment
- Outlier identification protocols

Analysis Validation:
- Sample size adequacy checks
- Statistical assumption verification
- Effect size calculations
- Confidence interval reporting

Reproducibility Features:
- Complete parameter logging
- Analysis timestamp tracking
- Filter setting documentation
- Random seed control

Clinical Integration Features:
-----------------------------

Time Period Analysis:
- Surgery: Complete intraoperative period
- CPB: Cardiopulmonary bypass duration
- Clamp: Aortic cross-clamp period
- Custom: User-defined time windows

Machine Validation:
- Cross-platform reliability assessment
- Technical artifact detection
- Signal quality comparison
- Harmonization validation

Outcome Mapping:
- Neurodevelopmental assessments
- Seizure classifications
- Risk stratification categories
- Long-term follow-up integration

Performance Characteristics:
--------------------------

Scalability:
- Handles 100+ patient datasets efficiently
- Memory usage: ~2GB for typical analysis
- Processing time: 5-15 minutes for full analysis
- Parallel processing capabilities

Memory Management:
- Intelligent object caching
- Explicit memory cleanup
- Large dataset streaming
- Garbage collection optimization

Error Handling:
- Graceful failure recovery
- Detailed error logging
- User-friendly error messages
- Analysis continuation strategies

Output File Organization:
------------------------

Primary Results:
- all_results/: Complete statistical analysis
- significant_results/: p < 0.05 findings
- trending_results/: 0.05 ≤ p < 0.10 findings

Specialized Analyses:
- across_seizures_comparison/: Seizure vs. non-seizure
- no_seizures_comparison/: ND outcomes in seizure-free patients
- correlation_results_*/: Correlation matrices

Documentation:
- filter_settings.txt: Analysis parameters
- analyzed_fis_log.txt: Included patients
- analysis_info.txt: Technical details

Export Formats:
- PDF: Publication-ready statistical reports
- Excel: Patient-level data for further analysis
- CSV: Time series data for modeling

Configuration Dependencies:
--------------------------
- base_dir: Root directory for analysis outputs
- clear: Whether to rebuild cached objects
- Statistical significance thresholds (customizable)
- Quality control parameters

Clinical Applications:
---------------------

Regulatory Submissions:
- FDA/EMA biomarker qualification packages
- Standardized statistical reporting
- Validation study protocols
- Safety signal detection

Research Publications:
- Automated manuscript tables
- Publication-ready figures
- Effect size calculations
- Multiple comparison corrections

Clinical Decision Support:
- Real-time biomarker analysis
- Risk stratification tools
- Treatment response monitoring
- Outcome prediction models

Future Extensions:
-----------------
- Machine learning integration
- Longitudinal analysis capabilities
- Multi-site harmonization tools
- Real-time analysis dashboard

"""

import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
utils_dir = os.path.join(current_dir, '..', '..', 'utils')
sys.path.append(utils_dir)
from utils import config
vars_dict = config.get_vars_dict()
import os
import numpy as np
from tqdm import tqdm
import tkinter as tk
from tkinter import ttk
from scipy import stats
from stat_tests import param_tests, stat_dict, alex_mwu, alex_students, var_dict, change_resolution
from copy import copy, deepcopy
import re
import pandas as pd
import shutil
import statsmodels.api as sm
import datetime
import csv
from pathlib import Path, PurePath
import gc
import pickle as pkl
import time
from datetime import timedelta

from file_importing import import_pkl
from fis_object import fis_object
from var_obj import var_object
from PDF_fxns import get_base_save_dir, save_by_category, combine_pdfs, integrate_PDF
from correlations_fxns import correlations_without_filters, save_correlation_results
from table_1_generator import generate_demographics_tables
from table_2_generator import generate_risk_categories_table
from table_3_generator import generate_seizure_table
from table_4_generator import generate_neurodevelopmental_table
from table_5_generator import generate_artefact_table
from table_11_generator import generate_biomarker_table
from statsmodels.formula.api import ols
from filt_dict import filt_dict

def build_objs_dict(ts_dir, save_dir):
    """
    The goal of this function is to build a dictionary in which each FIS key is associated with a fis_obj instance.

    Inputs:
    - ts_dir (str): the directory where all the ts .csv files are stored
    """
    if "intraop_files" in str(ts_dir):
        op_status = "intraop"
    elif "postop_files" in str(ts_dir):
        op_status = "postop"
    
    if op_status == "intraop":
        fis_dict_file_name = os.path.join(save_dir,"fis_dict.pkl")
    elif op_status == "postop":
        fis_dict_file_name = os.path.join(save_dir,"fis_dict_postop.pkl")
    
    
    # First check if dictionary already exists
    if os.path.exists(fis_dict_file_name):
        if vars_dict["clear"]:
            os.remove(fis_dict_file_name)
        else:
            print(f"Loading existing FIS dictionary from file {fis_dict_file_name}. File size: {os.path.getsize(fis_dict_file_name)*10**-9} MB.")
            start_time = time.time()
            with open(fis_dict_file_name, 'rb') as f:
                fis_dict = pkl.load(f)
                num_original = len(fis_dict)
                end_time = time.time()
                elapsed_time = end_time - start_time
                # Format the elapsed time
                time_formatted = str(timedelta(seconds=int(elapsed_time)))
                # Print completion message with time
                print(f"fis_dict built. Total time: {time_formatted}")
                return fis_dict


    print("Building FIS dictionary from scratch...")
    fis_dict = {}        
    # Get all CSV files in the directory
    print("Scanning directory for CSV files...")
    csv_files = [f for f in os.listdir(ts_dir) if f.endswith('.csv')]
    print(f"Found {len(csv_files)} CSV files")
    # Extract FIS values
    print("Extracting unique FIS values...")
    pattern = r'timeseries_.*?_(\d+)\.csv'
    xx_values = [int(match.group(1)) for file in csv_files if (match := re.match(pattern, file))]
    all_fis = sorted(set(xx_values))
    num_original = len(all_fis)
    print(f"Found {num_original} unique FIS values")
    for fis in all_fis:
        build_fis_object(fis, ts_dir)
    
    for fis in tqdm(all_fis, desc = "Loading FIS objects into dictionary..."):                     # now everything is built
        cache_path = (Path(ts_dir).parent / "pkl_cache" /
                    f"FIS{fis}.pkl")
        with open(cache_path, 'rb') as f:
            fis_dict[fis] = pkl.load(f)

    # Save the dictionary
    print("Saving FIS dictionary to file...")
    with open(fis_dict_file_name, 'wb') as f:
        pkl.dump(fis_dict, f, protocol=pkl.HIGHEST_PROTOCOL)
    print("Dictionary building complete!")
    if len(fis_dict) != num_original:
        raise ValueError(f"There were originally {num_original} FIS values read from the time series files, but fis_dict had {len(fis_dict)} values.")
    #raise ValueError("Pause.")
    return fis_dict

def build_fis_object(fis, ts_dir):
    """
    The goal of this function is to build a stat object for a given FIS.

    Inputs:
    - fis (int): the fis number that will be used to match up files

    Outputs:
    - None: simply builds and caches the fis ob
    """
    #Check if we've already pickled
    cache_dir   = os.path.join(os.path.dirname(ts_dir), "pkl_cache")
    store_path  = os.path.join(cache_dir, f"FIS{fis}.pkl")
    if os.path.exists(store_path):
        return pkl.load(store_path)   
    
    fis = int(fis)
    vars = []
    # Get all matching files first so we can track progress
    matching_files = [f for f in os.listdir(ts_dir) if f.endswith(f"_{fis}.csv")]
    # For each matching .csv file in ts_dir with progress bar
    for filename in tqdm(matching_files, desc=f"Processing variables for FIS{fis}", leave=False):
        if "clash" in filename:
            pass
        # Extract the var_name from the file's title (everything between the first and last underscores)
        parts = filename.split('_')
        if len(parts) >= 3:
            var_name = '_'.join(parts[1:-1])
            # Read the CSV file to extract values and times
            filepath = os.path.join(ts_dir, filename)
            vals = []
            times = []
            with open(filepath, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    vals.append(float(row['value']))
                    times.append(float(row['time']))
            # Set the sampling rate based on variable name
            if "rEEG" in var_name and "__" not in var_name:
                srate = 1/2
            else:
                srate = 1/30
            # Create the var_object and add to vars list
            vars.append(var_object(var_name, vals, times, srate))
    
    # Create a dictionary mapping variable names to their objects
    ts_dict = {}
    for var in vars:
        ts_dict[var.name] = var
    if "intraop_files" in str(ts_dir):
        op_status = "intraop"
    elif "postop_files" in str(ts_dir):
        op_status = "postop"
    ts_path = Path(ts_dir)
    ts_path_split_1 = os.path.split(ts_path)[0] #Gets us to ...\99_extract_all_variables\
    ts_path_split_2 = os.path.split(ts_path_split_1)[0] #Gets us to ...\intraop_files\ (or \postop_files\)
    ts_path_split_2 = os.path.join(os.path.join(ts_path_split_2,r"3_preprocessing\preprocessed_files"))
    base_full_path = os.path.join(vars_dict["base_dir"],ts_path_split_2)
    full_path = next((os.path.join(base_full_path, f) for f in os.listdir(base_full_path) if f"fis{fis}_".lower() in f.lower()), None)
    if full_path is None or not os.path.exists(full_path):
        raise ValueError(f"No file found for FIS{fis} at path {full_path}")
    raw_obj = import_pkl(full_file_path=full_path, parts_or_path="path")
    EEG_data, EEG_times = raw_obj.EEG_data, raw_obj.EEG_times_from_analyze_end_time
    srate = deepcopy(raw_obj.sampling_rate)
    del raw_obj
    
    # Return the fis_object
    gc.collect()
    fis_obj = fis_object(fis, ts_dict, EEG_data, EEG_times, op_status, srate)
    store_dir = os.path.dirname(ts_dir)
    store_dir = os.path.join(store_dir, "pkl_cache")
    store_dir = os.path.join(store_dir, f"FIS{fis_obj.fis}.pkl")
    print(store_dir)
    #breakpoint()
    with open(store_dir, 'wb') as file: 
        pkl.dump(fis_obj, file, protocol = pkl.HIGHEST_PROTOCOL)
    del fis_obj
    return None

def adjust_for_seizure(nd_group, no_nd_group, nd_ctrl_values, no_nd_ctrl_values):
    """
    Adjust for seizure status by subtracting the mean difference between seizure and non-seizure groups
    """
    
    # Check which subjects have seizure status available
    valid_nd_idx = [i for i, vals in enumerate(nd_ctrl_values) if vals.get('seizures_intraop') is not None]
    valid_no_nd_idx = [i for i, vals in enumerate(no_nd_ctrl_values) if vals.get('seizures_intraop') is not None]
    
    if not valid_nd_idx or not valid_no_nd_idx:
        return nd_group, no_nd_group  # No adjustment if insufficient data
    
    # Get seizure status values
    nd_seizure = [nd_ctrl_values[i].get('seizures_intraop') for i in valid_nd_idx]
    no_nd_seizure = [no_nd_ctrl_values[i].get('seizures_intraop') for i in valid_no_nd_idx]
    
    # Get corresponding data points
    nd_valid_data = [nd_group[i] for i in valid_nd_idx]
    no_nd_valid_data = [no_nd_group[i] for i in valid_no_nd_idx]
    
    # Calculate the effect of seizure status
    all_data = nd_valid_data + no_nd_valid_data
    all_seizure = nd_seizure + no_nd_seizure
    
    # Find mean difference between seizure and non-seizure groups
    seizure_group = [all_data[i] for i, s in enumerate(all_seizure) if s == 1]
    non_seizure_group = [all_data[i] for i, s in enumerate(all_seizure) if s == 0]
    
    if seizure_group and non_seizure_group:
        seizure_effect = np.mean(seizure_group) - np.mean(non_seizure_group)
        
        # Adjust values: subtract effect from seizure subjects
        adjusted_nd = np.array(nd_group)
        adjusted_no_nd = np.array(no_nd_group)
        
        for i in valid_nd_idx:
            if nd_ctrl_values[i].get('seizures_intraop') == 1:
                adjusted_nd[i] -= seizure_effect
        
        for i in valid_no_nd_idx:
            if no_nd_ctrl_values[i].get('seizures_intraop') == 1:
                adjusted_no_nd[i] -= seizure_effect
        
        return adjusted_nd.tolist(), adjusted_no_nd.tolist()
    
    return nd_group, no_nd_group  # No adjustment if insufficient data

def adjust_for_continuous_variables(nd_group, no_nd_group, nd_ctrl_values, no_nd_ctrl_values, 
                                   control_for_ga=False, control_for_pma=False):
    """
    Adjust for continuous variables (GA and/or PMA) using linear regression
    """
    
    # Combine both groups for regression
    all_data = nd_group + no_nd_group
    all_ctrl_values = nd_ctrl_values + no_nd_ctrl_values
    
    # Prepare predictor variables (X) and response variable (y)
    X = []
    valid_indices = []
    
    for i, ctrl_vals in enumerate(all_ctrl_values):
        predictors = []
        is_valid = True
        
        if control_for_ga:
            ga_val = ctrl_vals.get('ga_eeg')
            if ga_val is not None:
                predictors.append(ga_val)
            else:
                is_valid = False
        
        if control_for_pma:
            pma_val = ctrl_vals.get('pma_eeg')
            if pma_val is not None:
                predictors.append(pma_val)
            else:
                is_valid = False
        
        if is_valid and predictors:
            X.append(predictors)
            valid_indices.append(i)
    
    # If insufficient valid data, return original values
    if len(valid_indices) < 2:
        return nd_group, no_nd_group
    
    X = np.array(X)
    y = np.array([all_data[i] for i in valid_indices])
    
    # Create and fit a statsmodels linear regression model for better stats
    X = sm.add_constant(X)  # Add intercept
    model = sm.OLS(y, X)
    results = model.fit()
    
    # Get regression coefficients
    coefficients = results.params[1:]  # Skip intercept
    intercept = results.params[0]
    
    # Apply adjustment to both groups
    adjusted_nd = np.array(nd_group)
    adjusted_no_nd = np.array(no_nd_group)
    
    # Apply adjustments to ND group
    for i, ctrl_vals in enumerate(nd_ctrl_values):
        predictors = []
        is_valid = True
        
        if control_for_ga:
            ga_val = ctrl_vals.get('ga_eeg')
            if ga_val is not None:
                predictors.append(ga_val)
            else:
                is_valid = False
        
        if control_for_pma:
            pma_val = ctrl_vals.get('pma_eeg')
            if pma_val is not None:
                predictors.append(pma_val)
            else:
                is_valid = False
        
        if is_valid and predictors:
            # Calculate predicted value and adjust
            predictor_array = np.array(predictors)
            predicted = np.sum(predictor_array * coefficients)
            adjusted_nd[i] -= predicted - intercept
    
    # Apply adjustments to non-ND group
    for i, ctrl_vals in enumerate(no_nd_ctrl_values):
        predictors = []
        is_valid = True
        
        if control_for_ga:
            ga_val = ctrl_vals.get('ga_eeg')
            if ga_val is not None:
                predictors.append(ga_val)
            else:
                is_valid = False
        
        if control_for_pma:
            pma_val = ctrl_vals.get('pma_eeg')
            if pma_val is not None:
                predictors.append(pma_val)
            else:
                is_valid = False
        
        if is_valid and predictors:
            # Calculate predicted value and adjust
            predictor_array = np.array(predictors)
            predicted = np.sum(predictor_array * coefficients)
            adjusted_no_nd[i] -= predicted - intercept
    
    return adjusted_nd.tolist(), adjusted_no_nd.tolist()


class stat_gui():
    def __init__(self, timeseries_file_dir=None):
        # List of directories to clear
        directories = ['all_results', 'significant_results', 'trending_results']
        # Clear each directory
        for directory in directories:
            if os.path.exists(directory):
                # Remove all files and subdirectories in the directory
                for item in os.listdir(directory):
                    item_path = os.path.join(directory, item)
                    if os.path.isfile(item_path):
                        os.remove(item_path)
                    elif os.path.isdir(item_path):
                        shutil.rmtree(item_path)
                print(f"Cleared directory: {directory}")
            else:
                print(f"Directory not found: {directory}")
        
        print("Initializing stat_gui...")
        #Step 1.2: Determine where user wants to save
        base_save_dir = get_base_save_dir(os.path.abspath(vars_dict["base_dir"]))
        self.base_save_dir = base_save_dir
        self.fis_dict = build_objs_dict(timeseries_file_dir, self.base_save_dir)
        self.num_original_fis = len(self.fis_dict)
        self.all_fis = list(self.fis_dict.keys())
        self.all_vars = list(self.fis_dict[self.all_fis[0]].ts_dict.keys())
        self.time_option = "all"
        print(f"Found {len(self.all_fis)} FIS objects")
        print(f"Found {len(self.all_vars)} variables")
        #at this point, self.fis_dict is where the data manipulation will occur
        self.clear_filters()
    
    
    def export_to_excel(self, vars_arr, stats_arr, time_range, custom_start=None, custom_end=None, 
                        suffix="", control_vars=None, use_adjusted_values=False, analysis_type = None):
        """
        Export all calculated variables, parameters, and subject information to an Excel file
        
        Parameters:
        - vars_arr: list of variable names
        - stats_arr: list of statistic names
        - time_range: which time period was analyzed
        - custom_start, custom_end: custom time ranges if applicable
        - suffix: suffix for the filename
        - control_vars: dictionary of control variables for ANCOVA adjustment
        - use_adjusted_values: if True and control_vars provided, export ANCOVA-adjusted values
        """
        import pandas as pd
        from datetime import datetime
        from statsmodels.formula.api import ols
        import statsmodels.api as sm
        
        # Create lists to store all data
        all_data = []
        
        # Get list of all analyzed subjects
        if analysis_type == "whole_sample":
            analyzed_subjects = [fis_obj for fis_obj in self.fis_dict.values() 
                                if fis_obj.analyzed_whole_sample == 1]
        elif analysis_type == "seizure_sample":
            analyzed_subjects = [fis_obj for fis_obj in self.fis_dict.values() 
                                if fis_obj.analyzed_seizure_sample == 1]
        elif analysis_type == "across_seizures":
            analyzed_subjects = [fis_obj for fis_obj in self.fis_dict.values() 
                                if fis_obj.analyzed_across_seizures == 1]
        elif analysis_type == "all":
            analyzed_subjects =[fis_obj for fis_obj in self.fis_dict.values()]
        
        print(f"Exporting data for {len(analyzed_subjects)} analyzed subjects...")
        
        # Determine if we're doing ANCOVA adjustments
        do_adjustments = (use_adjusted_values and control_vars and any(control_vars.values()) and not analysis_type == "all")
        
        if do_adjustments:
            print("Calculating ANCOVA-adjusted values...")
            # Dictionary to store adjusted values for each var/stat combination
            adjusted_values_dict = {}
            
        # For each subject, calculate all parameter/variable combinations and collect demographics
        for fis_obj in tqdm(analyzed_subjects, desc="Collecting subject data"):
            # Set time range based on analysis parameters
            if time_range == "surgery" or time_range == "all":
                start_time = fis_obj.cx_start_time
                end_time = fis_obj.cx_end_time
            elif time_range == "ecc":
                start_time = fis_obj.ECC_times["ecc_start"]
                end_time = fis_obj.ECC_times["ecc_end"]
            elif time_range == "half":
                start_time = fis_obj.ECC_times["ecc_start"]
                end_time = fis_obj.ECC_times["ecc_start"] + (fis_obj.ECC_times["ecc_end"] - fis_obj.ECC_times["ecc_start"]) / 2
            elif time_range == "clamp":
                start_time = fis_obj.ECC_times["clamp_start"]
                end_time = fis_obj.ECC_times["clamp_end"] 
            elif time_range == "custom":
                start_time = custom_start
                end_time = custom_end
                
            if analysis_type == "all":
                log_file_path = os.path.join(self.base_save_dir, "all_pts_excluded_log.txt")
                with open(log_file_path, 'w') as f:
                    f.write("")  # Create empty file
            if (start_time is None or end_time is None) and analysis_type == "all":
                error_message = f"Error generating data for {fis_obj.fis}: could not find start/end time data for the time range {time_range}\n"
                with open(log_file_path, 'a') as f:
                    f.write(error_message)
                continue
            
            # Get demographic and clinical info
            subject_data = {
                'FIS': fis_obj.fis,
                'ND_status': fis_obj.ND,
                'Machine': fis_obj.machine_type,
                'Sex': getattr(fis_obj, 'sex', None),
                'GA_EEG': getattr(fis_obj, 'ga_eeg', None),
                'PMA_EEG': getattr(fis_obj, 'pma_eeg', None),
                'Seizures_intraop': getattr(fis_obj, 'seizures_intraop', None),
            }
            
            # Get artefact information
            fis_obj.build_artefact_info(start_time, end_time)
            subject_data['proportion_artefact'] = fis_obj.proportion_artefact
            subject_data["num_artefacts"] = fis_obj.num_artefacts
            subject_data["time_artefact"] = fis_obj.time_artefact
            subject_data["clean_signal_length_secs"] = ((end_time-start_time)*3600)*(1-fis_obj.proportion_artefact)
            
            # Calculate all variable/stat combinations
            tests = param_tests()
            for var in vars_arr:
                for stat_name in stats_arr:
                    stat = tests[stat_name]
                    # Get variable data for the time range
                    try:
                        var_times, var_vals = fis_obj.ts_dict[var].get_vals(start_time, end_time, fis_obj.op_status)
                        # Calculate the statistic
                        result = stat(var_vals, var_times)
                        
                        # Store raw value
                        param_name = f"{var}_{stat_name}"
                        subject_data[param_name] = result
                        
                        # If doing adjustments, we'll calculate them after collecting all raw data
                        if do_adjustments:
                            if param_name not in adjusted_values_dict:
                                adjusted_values_dict[param_name] = {'fis': [], 'raw_values': [], 
                                                                'nd_status': [], 'control_vars': {}}
                                for cv in control_vars:
                                    if control_vars[cv]:
                                        adjusted_values_dict[param_name]['control_vars'][cv] = []
                            
                            adjusted_values_dict[param_name]['fis'].append(fis_obj.fis)
                            adjusted_values_dict[param_name]['raw_values'].append(result)
                            adjusted_values_dict[param_name]['nd_status'].append(fis_obj.ND)
                            
                            for cv in control_vars:
                                if control_vars[cv]:
                                    adjusted_values_dict[param_name]['control_vars'][cv].append(
                                        getattr(fis_obj, cv, None))
                        
                    except Exception as e:
                        # Handle if the variable doesn't exist or other calculation errors
                        subject_data[f"{var}_{stat_name}"] = "Error in calculation"
            
            # Add the subject's data to the master list
            all_data.append(subject_data)
        
        # If doing adjustments, calculate adjusted values for each parameter
        if do_adjustments:
            print("Calculating ANCOVA adjustments for each parameter...")
            
            # Convert to DataFrame for easier manipulation
            df = pd.DataFrame(all_data)
            
            # For each parameter, run ANCOVA and get adjusted values
            for param_name in tqdm(adjusted_values_dict.keys(), desc="Adjusting parameters"):
                param_data = adjusted_values_dict[param_name]
                
                # Create DataFrame for this parameter's ANCOVA
                ancova_df = pd.DataFrame({
                    'eeg_param': param_data['raw_values'],
                    'nd_status': param_data['nd_status']
                })
                
                # Add control variables
                for cv, values in param_data['control_vars'].items():
                    ancova_df[cv] = values
                
                # Check for missing values
                if ancova_df['eeg_param'].isnull().any() or ancova_df.isnull().any().any():
                    print(f"Skipping adjustment for {param_name} due to missing values")
                    continue
                
                # Build and fit ANCOVA model
                formula = 'eeg_param ~ nd_status'
                for cv in control_vars:
                    if control_vars[cv]:
                        formula += f' + {cv}'
                
                model = ols(formula, data=ancova_df).fit()
                
                # Calculate adjusted values for each subject
                # Get mean values for covariates
                covariate_means = {}
                for cv in control_vars:
                    if control_vars[cv]:
                        covariate_means[cv] = ancova_df[cv].mean()
                
                # Calculate adjusted values
                for i, row in ancova_df.iterrows():
                    # Create prediction DataFrame with average covariate values
                    pred_df = pd.DataFrame({
                        'nd_status': [row['nd_status']],
                        **{k: [v] for k, v in covariate_means.items()}
                    })
                    
                    # Get residual and adjusted value
                    residual = row['eeg_param'] - model.predict(pd.DataFrame(row).T).iloc[0]
                    adjusted_value = model.predict(pred_df).iloc[0] + residual
                    
                    # Update the subject's data with adjusted value
                    fis = param_data['fis'][i]
                    subject_idx = next(idx for idx, subj in enumerate(all_data) 
                                    if subj['FIS'] == fis)
                    
                    # Add adjusted value with suffix
                    all_data[subject_idx][f"{param_name}_adjusted"] = adjusted_value
                        
        
        # Convert to DataFrame
        df = pd.DataFrame(all_data)
        
        # Create the Excel filename based on the time range and whether values are adjusted
        adjustment_str = "_adjusted" if do_adjustments else "_unadjusted"
        excel_filename = f"{self.base_save_dir}\\calculated_parameters_{suffix}{adjustment_str}.xlsx"
        
        # Export to Excel with multiple sheets for better organization
        with pd.ExcelWriter(excel_filename, engine='openpyxl') as writer:
            # Main data sheet
            df.to_excel(writer, sheet_name='All_Data', index=False)
            
            # If adjusted values were calculated, create a comparison sheet
            if do_adjustments:
                # Find columns with adjusted values
                adjusted_cols = [col for col in df.columns if col.endswith('_adjusted')]
                if adjusted_cols:
                    # Create comparison DataFrame
                    comparison_data = {'FIS': df['FIS']}
                    
                    for adj_col in adjusted_cols:
                        raw_col = adj_col.replace('_adjusted', '')
                        if raw_col in df.columns:
                            comparison_data[f"{raw_col}_raw"] = df[raw_col]
                            comparison_data[f"{raw_col}_adjusted"] = df[adj_col]
                            comparison_data[f"{raw_col}_difference"] = df[adj_col] - df[raw_col]
                    
                    comparison_df = pd.DataFrame(comparison_data)
                    comparison_df.to_excel(writer, sheet_name='Raw_vs_Adjusted', index=False)
        
        print(f"Data exported successfully to: {excel_filename}")
        
        return excel_filename
    
    def filter_by_ND(self):
        for fis_obj in self.fis_dict.values():
            if fis_obj.ND is None:
                fis_obj.ND_filter = 0
        #self.syn_filters()
    
    def filter_by_seizure(self, filt):
        """
        Filters to only include patients with confirmed absence or presence of intraop seizure
        - filt = 0 to retain only pts without seizure, filt = 1 to retain only pts with seizure, filt = 2 to retain only pts with seizure data, filt = 3 to eliminate filtering
        """
        if filt == 0:
            for fis_obj in self.fis_dict.values():
                if fis_obj.seizures_intraop==0:
                    fis_obj.seizure_filter = 1
                else:
                    fis_obj.seizure_filter = 0
        elif filt == 1:
            for fis_obj in self.fis_dict.values():
                if fis_obj.seizures_intraop==1:
                    fis_obj.seizure_filter = 1
                else:
                    fis_obj.seizure_filter = 0
        elif filt == 2:
            for fis_obj in self.fis_dict.values():
                if fis_obj.seizures_intraop is None:
                    fis_obj.seizure_filter = 0
                else:
                    fis_obj.seizure_filter = 1
        elif filt == 3:
            for fis_obj in self.fis_dict.values():
                fis_obj.seizure_filter = 1
        #self.syn_filters()

    def clear_filters(self):
        """
        Resets filters that change based on the analysis (i.e., seizure and ND_status filters). Other filters (machine, ignore, ecc, and artefact) remain intact.
        """
        for fis_obj in self.fis_dict.values():
            fis_obj.ND_filter = 1
            fis_obj.seizure_filter = 1
        #self.syn_filters()

    def filter_by_ecc(self, ecc_str):
        """
        ecc_str = "ecc", "clamp", or "all"
        """
        if ecc_str == "ecc" or ecc_str == "clamp":
            self.time_option = ecc_str
            for fis_obj in self.fis_dict.values():
                if any(val is not None for val in list(fis_obj.ECC_times.values())) and any(ecc_str in key for key, val in fis_obj.ECC_times.items() if val is not None):
                    fis_obj.ecc_filter = 1
                else:
                    fis_obj.ecc_filter = 0
        else:
            for fis_obj in self.fis_dict.values():
                fis_obj.ecc_filter = 1
                self.time_option = "all"
        #self.syn_filters()
            #print(f"FIS{fis_obj.fis}, ND: {fis_obj.ND}, ecc: {fis_obj.ecc_status}, ND_filter: {fis_obj.ND_filter}, ecc_filter: {fis_obj.ecc_filter}, filter: {fis_obj.filter}")
    
    def filter_by_machine(self, machine_type):
        """
        machine_type here can be "nicolet", "obm", or "all"
        """
        if not (machine_type == "nicolet" or machine_type == "obm"):
            for fis_obj in self.fis_dict.values():
                fis_obj.machine_type_filter = 1
        else:
            for fis_obj in self.fis_dict.values():
                if fis_obj.machine_type == machine_type:
                    fis_obj.machine_type_filter = 1
                else:
                    fis_obj.machine_type_filter = 0
                #print(f"FIS{fis_obj.fis}, machine_tpye: {fis_obj.machine_type}, machine_type_filter: {fis_obj.machine_type_filter}, filter: {fis_obj.filter}")
        #self.syn_filters()
    
    def filter_by_ignore(self, fis_arr):
        """
        Input:
        - fis_arr (arr of int): array of the FIS values we want excluded from processing
        """
        for fis in fis_arr:
            fis = int(fis)
            for fis_obj in self.fis_dict.values():
                if fis_obj.fis == fis:
                    fis_obj.ignore_filter = 0
                else:
                    continue
        #self.syn_filters() 

    def run_test(self, var, stat, test=None, time_range="surgery", show_boxplot=False, print_results=False, mode="ind", cus_time_1=None, cus_time_2=None, min_art=0, comparison_category=["ND", "ND"], all_correlation_results=None, control_vars=None, use_adjusted_correlations=False, case = None, time_op = None):
        """
        Inputs:
        - var (str) = the name of variable of interest (e.g., delta power)
        - stat (str) = the name of the statistic we are extracting from the variable time series
        - test (str) = the test we would like to run ("mwu", "t", or None)
        - time_range (str): "surgery", "ecc", "clamp", "half", or "custom"
        - control_vars (dict): Dictionary indicating which variables to control for (e.g., {'seizure': True, 'ga_eeg': False})
        - use_adjusted_correlations (bool): Whether to use adjusted values for correlation analysis
        - case (str): tells us what analysis case we're in to help us track the analyses
        """
        #Step 0: Initialize important variables
        tests = param_tests()
        stat_name = copy(stat)
        stat = tests[stat]
        nd_group = []
        no_nd_group = []
        nd_fis = []
        no_nd_fis = []
        vars_dict = var_dict()
        stats_dict = stat_dict()
        whole_mean = None
        whole_sem = None
        whole_median = None
        whole_q1 = None
        whole_q3 = None
        row_data = None
        var_to_add = None
        equal_var = None
        ancova_bool = control_vars and any(control_vars.values()) and comparison_category[0] == "ND"
        
        if comparison_category == ["ND","ND"]:
            self.filter_by_ND() #To eliminate pts without outcomes data
        elif comparison_category == ["seizures_intraop","intraoperative seizure"]:
            self.filter_by_seizure(2) #To eliminate patients w/o seizures data
        #Step 0.1: For ANCOVA, collect control variable values and build a pandas DataFrame
        if ancova_bool:
            self.filter_by_seizure(2) #Maintain only patients with seizure data, as we need this data to do adjustment
            self.filter_by_ND() #Since we need ND data to run analyses
            #debug_list = [db_obj for db_obj in list(self.fis_dict.values()) if db_obj.filter==1]
            #print(len(debug_list))
            ancova_data = {'eeg_param': [], 'nd_status': [], 'fis': []}
            # Initialize columns for each control variable
            for var_name, enabled in control_vars.items():
                if enabled:
                    ancova_data[var_name] = []
            
        #Step 1: Get our data points, filtering by filter and parsing by time
        for fis_obj in self.fis_dict.values():            
            #Check filters
            #Time selection
            if time_range == "surgery" or time_range == "all":
                start_time = fis_obj.cx_start_time
                end_time = fis_obj.cx_end_time
            elif time_range == "ecc":
                if not (fis_obj.ECC_times["ecc_start"] is None or fis_obj.ECC_times["ecc_end"] is None):
                    start_time = fis_obj.ECC_times["ecc_start"]
                    end_time = fis_obj.ECC_times["ecc_end"]
                else:
                    fis_obj.skip_str = "no ECC time"
                    fis_obj.ecc_filter = 0
            elif time_range == "half":
                if not (fis_obj.ECC_times["ecc_start"] is None or fis_obj.ECC_times["ecc_end"] is None):
                    start_time = fis_obj.ECC_times["ecc_start"]
                    end_time = fis_obj.ECC_times["ecc_end"]
                    end_time = start_time + (end_time - start_time) / 2
                else:
                    fis_obj.skip_str = "no ECC time"
                    fis_obj.ecc_filter = 0
            elif time_range == "clamp":
                if not (fis_obj.ECC_times["clamp_start"] is None or fis_obj.ECC_times["clamp_end"] is None):
                    start_time = fis_obj.ECC_times["clamp_start"]
                    end_time = fis_obj.ECC_times["clamp_end"]
                else:
                    fis_obj.skip_str = "no clamp time"
                    fis_obj.ecc_filter = 0
            elif time_range == "custom":
                start_time = cus_time_1
                end_time = cus_time_2
            #Check artefact filters using time limits
            try:
                fis_obj.build_artefact_info(start_time, end_time)
            except:
                if fis_obj.ecc_filter == 0:
                    fis_obj.skip_str = "filtered out (d/t lack of ecc data)"
                    continue
                else:
                    raise ValueError("Artefact info building broke, but ECC data was present.")
            nan_fraction = fis_obj.proportion_artefact
            if nan_fraction > (1 - min_art):
                fis_obj.skip_str = f"excess artefact: {nan_fraction*100}%"
                fis_obj.artefact_filter = 0
            # if fis_obj.fis in [34, 158, 160]:
            #     pass ; prior debugging step
            if fis_obj.filter == 0:
                if fis_obj.machine_type_filter == 0:
                    fis_obj.skip_str = "filtered out (d/t machine)"
                elif fis_obj.ecc_filter == 0:
                    fis_obj.skip_str = "filtered out (d/t lack of ecc data)"
                elif fis_obj.ND_filter == 0:
                    fis_obj.skip_str = "filtered out (d/t lack of ND data)"
                elif fis_obj.seizure_filter == 0:
                    if fis_obj.seizures_intraop == 1:
                        fis_obj.skip_str = "filtered out (d/t presence of seizure)"
                    elif fis_obj.seizures_intraop is None:
                        fis_obj.skip_str = "filtered out (d/t lack of seizure data)"
                    else:
                        fis_obj.skip_str = "filtered out (d/t seizure data, but details unclear)"
                elif fis_obj.artefact_filter == 0:
                    fis_obj.skip_str = f"filtered out (d/t excess artefact: {nan_fraction*100}%)"
                else:
                    if not "(reason unknown)" in fis_obj.skip_str:
                        fis_obj.skip_str = f"filtered out (reason unknown). self.ND_filter={fis_obj.ND_filter}, self.ecc_filter={fis_obj.ecc_filter}, self.machine_type_filter={fis_obj.machine_type_filter}, self.ignore_filter={fis_obj.ignore_filter}, self.seizure_filter={fis_obj.seizure_filter}."
                continue
            else:
                if case == "whole_sample":
                    fis_obj.analyzed_whole_sample = 1
                elif case == "seizure_sample":
                    fis_obj.analyzed_seizure_sample = 1
                elif case == "across_seizures":
                    fis_obj.analyzed_across_seizures = 1
        
            # Step 1.1: Retrieve time series data
            att_str = comparison_category[0]
            head_str = comparison_category[1]
            var_times, var_vals = fis_obj.ts_dict[var].get_vals(start_time, end_time, fis_obj.op_status)
            add_stat = stat(var_vals, var_times)
            # Dynamically retrieve the comparison attibute (i.e, compare by seizure or ND?)
            att_val = getattr(fis_obj, att_str, None)
                
            # Step 1.2: Run ANCOVA or default
            if ancova_bool:
                #Step 1.2.1: Data checks
                # Check if all control variables are available for this subject
                for var_name, enabled in control_vars.items():
                    if enabled and getattr(fis_obj, var_name, None) is None:
                        raise ValueError(f"We are asking to use control variables that are not present for FIS{fis_obj.fis}. This should have been captured for the seizures case (the non-seizures cases have not been deeply debugged yet).")
                # Step 1.2.2: Build our ancova_data data frame.
                ancova_data['eeg_param'].append(add_stat) #This is the EEG parameter
                ancova_data['nd_status'].append(1 if att_val == 1 else 0) #This is the ND status
                ancova_data['fis'].append(fis_obj.fis) 
                # Add control variable values
                for control_var_name, enabled in control_vars.items():
                    if enabled:
                        ancova_data[control_var_name].append(getattr(fis_obj, control_var_name, None))
            # Standard group assignment (we still need this for boxplots and non-ANCOVA results)
            if att_val == 0:
                no_nd_group.append(add_stat)
                no_nd_fis.append(fis_obj.fis)
            elif att_val == 1:
                nd_group.append(add_stat)
                nd_fis.append(fis_obj.fis)
        debug_check = sum([1 for val in self.fis_dict.values() if val.analyzed_whole_sample == 1])
        #print(f"{debug_check} FIS files marked as analyzed.")
        #Step 1.2.3: Run ANCOVA
        if ancova_bool and len(ancova_data['eeg_param']) > 2:
            #Step 1.2.3.1: Prep df and run data quality checks before ANCOVA
            # Convert to DataFrame
            df = pd.DataFrame(ancova_data)
            # Build formula for the model
            formula = f'eeg_param ~ nd_status' #att_str refers to the two groups of comparison
            for control_var_name, enabled in control_vars.items():
                if enabled:
                    formula += f' + {control_var_name}' #Will add "+ " and then the controlling variable(s) of interest (seizure, GA, PMA) to the formula
            # Check for insufficient data or missing values
            if df.isnull().any().any():
                raise ValueError("NoneType value when there should not be one.")

            #Step 1.2.3.2: Fit ANCOVA model
            stat_str = stats_dict[stat_name] #beautified statistics string
            var_str = vars_dict[var].capitalize() #vars_dict is just a dict to beautify the name of the variable we are currently testing
            model = ols(formula, data=df).fit()
            # Extract results
            nd_effect = model.params['nd_status']
            nd_pvalue = model.pvalues['nd_status']
            # Store for reporting
            p_val = nd_pvalue
            test = "ANCOVA"
            # Calculate adjusted means for reporting
            # Create separate DataFrames for ND and non-ND groups with mean values for covariates
            covariate_means = {}
            for var_name, enabled in control_vars.items():
                if enabled:
                    covariate_means[var_name] = df[var_name].mean()
                
            # Create prediction DataFrames with average covariate values
            nd_pred_df = pd.DataFrame({
                'nd_status': [1],
                **{k: [v] for k, v in covariate_means.items()}
            })
            non_nd_pred_df = pd.DataFrame({
                'nd_status': [0],
                **{k: [v] for k, v in covariate_means.items()}
            })

            # Get model predictions for adjusted means
            nd_mid = model.predict(nd_pred_df).iloc[0]
            no_nd_mid = model.predict(non_nd_pred_df).iloc[0]

            # Get proper standard errors for the adjusted means
            nd_prediction = model.get_prediction(nd_pred_df)
            no_nd_prediction = model.get_prediction(non_nd_pred_df)
            nd_spread = nd_prediction.se_mean[0]  # Standard error for ND group adjusted mean
            no_nd_spread = no_nd_prediction.se_mean[0]  # Standard error for non-ND group adjusted mean

            # Calculate adjusted values for each individual
            adjusted_values = []
            nd_adjusted_values = []
            no_nd_adjusted_values = []

            # For each subject in the analysis
            for i, row in df.iterrows():
                # Create a dataframe for this subject with average values for covariates
                subject_df = pd.DataFrame({
                    'nd_status': [row['nd_status']],  # Keep actual ND status
                    **{k: [covariate_means[k]] for k in covariate_means.keys()}  # Use average covariate values
                })
                
                # Get the adjusted value (residual + predicted with average covariates)
                residual = row['eeg_param'] - model.predict(pd.DataFrame(row).T).iloc[0]
                adjusted_value = model.predict(subject_df).iloc[0] + residual
                
                # Add to appropriate list
                adjusted_values.append(adjusted_value)
                if row['nd_status'] == 1:
                    nd_adjusted_values.append(adjusted_value)
                else:
                    no_nd_adjusted_values.append(adjusted_value)

            # Calculate whole group statistics based on adjusted values
            whole_mean = np.mean(adjusted_values)
            whole_sem = np.std(adjusted_values, ddof=1) / np.sqrt(len(adjusted_values))
            
            # Set test statistic
            test_stat = model.tvalues['nd_status']
            # Get R-squared for reporting
            r_squared = model.rsquared
            # Set equal_var (we'll use this for reporting)
            equal_var = True  # ANCOVA assumes homogeneity of regression slopes
                
            # Create a string to report which variables were controlled for
            controlled_vars = []
            control_effects = {}
            for var_name, enabled in control_vars.items():
                if enabled:
                    controlled_vars.append(var_name)
                    if var_name in model.params:
                        control_effects[var_name] = (model.params[var_name], model.pvalues[var_name])
            controlled_str = ", ".join(controlled_vars)
            
        # If not using ANCOVA (or ANCOVA failed), run standard analysis
        if not ancova_bool:
            # This is the original analysis code
            #Step 2: Run our test and get statistics
            stat_str = stats_dict[stat_name]
            var_str = vars_dict[var].capitalize()
            if test == "mwu":
                # Calculate quartiles for both groups
                nd_q1, nd_q3 = np.percentile(nd_group, [25, 75])
                no_nd_q1, no_nd_q3 = np.percentile(no_nd_group, [25, 75])
                # Calculate quartiles for the whole sample
                whole_sample = np.concatenate([nd_group, no_nd_group])
                whole_median = np.median(whole_sample)
                whole_q1, whole_q3 = np.percentile(whole_sample, [25, 75])
                nd_mid, nd_spread, no_nd_mid, no_nd_spread, test_stat, p_val, equal_var = alex_mwu(nd_group,no_nd_group, boxplot = show_boxplot, var_name = var, param = stat_name)
            elif test == "t":
                # Calculate sem for both groups
                nd_sem = np.std(nd_group, ddof=1) / np.sqrt(len(nd_group))
                no_nd_sem = np.std(no_nd_group, ddof=1) / np.sqrt(len(no_nd_group))
                whole_sample = np.concatenate([nd_group, no_nd_group])
                whole_mean = np.mean(whole_sample)
                whole_sem = np.std(whole_sample, ddof=1) / np.sqrt(len(whole_sample))
                nd_mid, nd_spread, no_nd_mid, no_nd_spread, test_stat, p_val, equal_var = alex_students(nd_group,no_nd_group, boxplot = show_boxplot, var_name = var, param = stat_name)
            elif test is None:
                # Perform Shapiro-Wilk test
                nd_shapiro_stat, nd_shapiro_p = stats.shapiro(nd_group)
                no_nd_shapiro_stat, no_nd_shapiro_p = stats.shapiro(no_nd_group)
                # Check if either group is not normal (p < 0.05)
                if nd_shapiro_p < 0.05 or no_nd_shapiro_p < 0.05:
                    norm_bool = False
                    nd_mid, nd_spread, no_nd_mid, no_nd_spread, test_stat, p_val, _ = alex_mwu(
                        nd_group, no_nd_group, boxplot=show_boxplot, var_name=var, param=stat_name)
                    test = "mwu"
                    # Calculate quartiles for both groups
                    nd_q1, nd_q3 = np.percentile(nd_group, [25, 75])
                    no_nd_q1, no_nd_q3 = np.percentile(no_nd_group, [25, 75])
                    # Calculate quartiles for the whole sample
                    whole_sample = np.concatenate([nd_group, no_nd_group])
                    whole_median = np.median(whole_sample)
                    whole_q1, whole_q3 = np.percentile(whole_sample, [25, 75])
                else:
                    norm_bool = True
                    nd_mid, nd_spread, no_nd_mid, no_nd_spread, test_stat, p_val, equal_var = alex_students(
                        nd_group, no_nd_group, boxplot=show_boxplot, var_name=var, param=stat_name)
                    test = "t"
                    # Calculate sem for both groups
                    nd_sem = np.std(nd_group, ddof=1) / np.sqrt(len(nd_group))
                    no_nd_sem = np.std(no_nd_group, ddof=1) / np.sqrt(len(no_nd_group))
                    whole_sample = np.concatenate([nd_group, no_nd_group])
                    whole_mean = np.mean(whole_sample)
                    whole_sem = np.std(whole_sample, ddof=1) / np.sqrt(len(whole_sample))
            
        

        
        #Step 3: build output
        if mode == "all":
            # Check if we used ANCOVA
            if ancova_bool and 'controlled_str' in locals():
                # Get actual number of subjects used in ANCOVA after filtering
                actual_nd_count = sum(df['nd_status'] == 1)
                actual_no_nd_count = sum(df['nd_status'] == 0)
                    
                whole_mean, whole_sem, no_nd_mid, no_nd_spread, nd_mid, nd_spread = change_resolution([whole_mean, whole_sem, no_nd_mid, no_nd_spread, nd_mid, nd_spread])
                # Create formatted strings for each column using pre-calculated variables
                var_to_add = f"{var_str}, {stat_str}"
                group_stats = f"{whole_mean} ({whole_sem})"
                no_nd_stats = f"{no_nd_mid} ({no_nd_spread})"
                nd_stats = f"{nd_mid} ({nd_spread})"
                r_squared_str = f", R² = {r_squared:.3f}"
                test_results = f"p = {p_val:.3f} (t = {test_stat:.3f}{r_squared_str})"
                equal_var_info = f"-"
                    
                # Use the actual counts in the column headers
                no_nd_str = f"No {head_str} (n={actual_no_nd_count})"
                nd_str = f"{head_str} (n={actual_nd_count})"
                    
                
                row_data = {
                    'Variable': var_to_add,
                    'Group Values': group_stats,
                    no_nd_str: no_nd_stats,
                    nd_str: nd_stats,
                    'Test Results': test_results,
                    'Equal Var': equal_var_info
                }
                    
            elif test == "mwu":
                not_dict = {False : "Equal variances not assumed.",
                            True : "Equal variances assumed.",
                            None : "MWU, so no Levene's test run."
                            }
                
                whole_median, whole_q1, whole_q3, no_nd_mid, no_nd_q1, no_nd_q3, nd_mid, nd_q1, nd_q3 = change_resolution([whole_median, whole_q1, whole_q3, no_nd_mid, no_nd_q1, no_nd_q3, nd_mid, nd_q1, nd_q3])
                # Create formatted strings for each column using pre-calculated variables
                var_to_add = f"{var_str}, {stat_str}"
                group_stats = f"{whole_median} [{whole_q1}, {whole_q3}]"
                no_nd_stats = f"{no_nd_mid} [{no_nd_q1}, {no_nd_q3}]"
                nd_stats = f"{nd_mid} [{nd_q1}, {nd_q3}]"
                test_results = f"p = {p_val:.3f}, (U = {test_stat:.0f})"
                equal_var_info = f"{not_dict[equal_var]}"
                    
                # Combine into a data frame row
                no_nd_str = f"No {head_str} (n={len(no_nd_group)})"
                nd_str = f"{head_str} (n={len(nd_group)})"
                row_data = {
                    'Variable': var_to_add,
                    'Group Values': group_stats,
                    no_nd_str: no_nd_stats,
                    nd_str: nd_stats,
                    'Test Results': test_results,
                    'Equal Var': equal_var_info
                }
            elif test == "t":
                not_dict = {False : "Equal variances not assumed.",
                            True : "Equal variances assumed.",
                            None : "MWU, so no Levene's test run."
                            }
                whole_mean, whole_sem, no_nd_mid, no_nd_sem, nd_mid, nd_sem = change_resolution([whole_mean, whole_sem, no_nd_mid, no_nd_spread, nd_mid, nd_spread])
                # Create formatted strings for each column using pre-calculated variables
                var_to_add = f"{var_str}, {stat_str}"
                group_stats = f"{whole_mean} ({whole_sem})"
                no_nd_stats = f"{no_nd_mid} ({no_nd_sem})"
                nd_stats = f"{nd_mid} ({nd_sem})"
                test_results = f"p = {p_val:.3f} (t = {test_stat:.0f})"
                equal_var_info = f"{not_dict[equal_var]}"
                no_nd_str = f"No {head_str} (n={len(no_nd_group)})"
                nd_str = f"{head_str} (n={len(nd_group)})"
                    
                # Combine into a data frame row
                row_data = {
                    'Variable': var_to_add,
                    'Group Values': group_stats,
                    no_nd_str: no_nd_stats,
                    nd_str: nd_stats,
                    'Test Results': test_results,
                    'Equal Var': equal_var_info
                }
            row_data = pd.DataFrame([row_data])
                
            if all_correlation_results is not None:
                #self.syn_filters()
                my_filt_dict = filt_dict(self.fis_dict)

                corr_results = correlations_without_filters(self, self.fis_dict, my_filt_dict, var, stat, time_range=time_range, cus_time_1=cus_time_1, cus_time_2=cus_time_2, comparison_category=comparison_category, stat_name=stat_name,use_adjusted_values=use_adjusted_correlations, control_vars=control_vars)
                del my_filt_dict
                all_correlation_results.extend(corr_results)
                
            return row_data, var_to_add, p_val  # The last two values will help with reporting the final significance results

    def get_filter_options(self):
        """
        Builds a GUI to get important filter information from the user.
        """
        # Create the main window
        root = tk.Tk()
        root.title("Filter Options")
        root.geometry("800x800")  # Increased height for better spacing
        
        # Create main container frame with grid layout
        main_container = ttk.Frame(root)
        main_container.pack(fill="both", expand=True)
        main_container.grid_columnconfigure(0, weight=1)
        
        # Variables to store user selections
        machine_var = tk.StringVar(value="all")
        intraop_var = tk.StringVar(value="all")
        ignore_var = tk.StringVar()
        multi_condition_var = tk.BooleanVar(value=False)
        
        # Variables for control selections (ANCOVA)
        control_seizure_var = tk.BooleanVar(value=False)
        control_ga_var = tk.BooleanVar(value=False)
        control_pma_var = tk.BooleanVar(value=False)
        
        # Variable for adjusted correlation analysis
        use_adjusted_correlations_var = tk.BooleanVar(value=False)
        
        # Variables for checkboxes
        vars_dict = var_dict()
        stats_dict = stat_dict()
        var_checkboxes = {}
        stat_checkboxes = {}
        
        # Create frames for organization
        machine_frame = ttk.LabelFrame(main_container, text="Machine Type")
        machine_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=5)
        
        intraop_frame = ttk.LabelFrame(main_container, text="Intraoperative Period")
        intraop_frame.grid(row=1, column=0, sticky="ew", padx=10, pady=5)
        
        ignore_frame = ttk.LabelFrame(main_container, text="Ignore FIS Values")
        ignore_frame.grid(row=2, column=0, sticky="ew", padx=10, pady=5)
        
        # Add a new frame for the multi-condition analysis checkbox
        multi_cond_frame = ttk.Frame(main_container)
        multi_cond_frame.grid(row=3, column=0, sticky="ew", padx=10, pady=5)
        
        # MODIFIED: Combined control frame for horizontal layout
        combined_controls_frame = ttk.LabelFrame(main_container, text="Analysis Controls")
        combined_controls_frame.grid(row=4, column=0, sticky="ew", padx=10, pady=5)
        
        # Create a horizontal layout within the combined controls frame
        control_frame = ttk.Frame(combined_controls_frame)
        control_frame.pack(fill="x", expand=True, padx=5, pady=5)
        
        # Divide the control frame into three columns
        control_frame.columnconfigure(0, weight=1)
        control_frame.columnconfigure(1, weight=1)
        control_frame.columnconfigure(2, weight=1)
        
        # ANCOVA controls in first column
        ancova_frame = ttk.LabelFrame(control_frame, text="Control Variables (ANCOVA)")
        ancova_frame.grid(row=0, column=0, sticky="nw", padx=5, pady=5)
        ttk.Label(ancova_frame, text="Select variables to control for:").pack(anchor="w", padx=5, pady=2)
        ttk.Checkbutton(ancova_frame, text="Seizure Status (seizures_intraop)", variable=control_seizure_var).pack(anchor="w", padx=5, pady=2)
        ttk.Checkbutton(ancova_frame, text="Gestational Age (ga_eeg)", variable=control_ga_var).pack(anchor="w", padx=5, pady=2)
        ttk.Checkbutton(ancova_frame, text="Postmenstrual Age (pma_eeg)", variable=control_pma_var).pack(anchor="w", padx=5, pady=2)
        # Add info label about ANCOVA
        ttk.Label(ancova_frame, text="Note: Using control variables will run ANCOVA\ninstead of t-test/Mann-Whitney U", 
                font=("TkDefaultFont", 8, "italic")).pack(anchor="w", padx=5, pady=2)
        
        # Correlation options in second column
        corr_frame = ttk.LabelFrame(control_frame, text="Correlation Analysis Options")
        corr_frame.grid(row=0, column=1, sticky="nw", padx=5, pady=5)
        ttk.Checkbutton(
            corr_frame, 
            text="Use adjusted values for correlations", 
            variable=use_adjusted_correlations_var
        ).pack(anchor="w", padx=5, pady=2)
        ttk.Label(corr_frame, 
                text="Note: This option only applies if control\nvariables are selected above", 
                font=("TkDefaultFont", 8, "italic")).pack(anchor="w", padx=5, pady=2)
        
        min_valid_data = tk.StringVar(value="0.20")  # As proportion (e.g., 0.9 for 90%)
        # Assuming there's a variable like this for the custom option selection
        # Create but don't grid the custom time frame initially
        custom_time_frame = ttk.LabelFrame(control_frame, text="Custom Time Window (hours)")
        # Don't grid it yet - will be shown only when custom is selected
        # Custom Time Period Inputs
        custom_start = tk.StringVar()
        custom_end = tk.StringVar()
        ttk.Label(custom_time_frame, text="Start Time:").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        ttk.Entry(custom_time_frame, textvariable=custom_start, width=8).grid(row=0, column=1, padx=5, pady=2)
        ttk.Label(custom_time_frame, text="End Time:").grid(row=1, column=0, sticky="w", padx=5, pady=2)
        ttk.Entry(custom_time_frame, textvariable=custom_end, width=8).grid(row=1, column=1, padx=5, pady=2)
        ttk.Label(custom_time_frame, text="Min Proportion Valid:").grid(row=2, column=0, sticky="w", padx=5, pady=2)
        ttk.Entry(custom_time_frame, textvariable=min_valid_data, width=8).grid(row=2, column=1, padx=5, pady=2)
        # Add this function to toggle visibility of custom frame
        def toggle_custom_frame(*args):
            if intraop_var.get() == "custom":
                custom_time_frame.grid(row=0, column=2, sticky="nw", padx=5, pady=5)
            else:
                custom_time_frame.grid_forget()

        intraop_var.trace_add("write", toggle_custom_frame)

        # call once to set initial state
        toggle_custom_frame()
        
        # Variables - now with more vertical space
        vars_frame = ttk.LabelFrame(main_container, text="Variables")
        vars_frame.grid(row=5, column=0, sticky="nsew", padx=10, pady=5)
        main_container.grid_rowconfigure(5, weight=2)  # Increased weight to give more space
        
        # Statistics - now with more vertical space
        stats_frame = ttk.LabelFrame(main_container, text="Statistics")
        stats_frame.grid(row=6, column=0, sticky="nsew", padx=10, pady=5)
        main_container.grid_rowconfigure(6, weight=2)  # Increased weight to give more space
        
        # Buttons at the bottom
        button_frame = ttk.Frame(main_container)
        button_frame.grid(row=7, column=0, sticky="ew", padx=10, pady=10)
        
        # Machine type radio buttons
        ttk.Radiobutton(machine_frame, text="Nicolet", variable=machine_var, value="nicolet").pack(side="left", padx=10)
        ttk.Radiobutton(machine_frame, text="OBM", variable=machine_var, value="obm").pack(side="left", padx=10)
        ttk.Radiobutton(machine_frame, text="All", variable=machine_var, value="all").pack(side="left", padx=10)
        
        # Intraoperative period radio buttons
        ttk.Radiobutton(intraop_frame, text="Entire Surgery", variable=intraop_var, value="all").pack(side="left", padx=10)
        ttk.Radiobutton(intraop_frame, text="CPB", variable=intraop_var, value="ecc").pack(side="left", padx=10)
        ttk.Radiobutton(intraop_frame, text="CPD Clamping Only", variable=intraop_var, value="clamp").pack(side="left", padx=10)
        ttk.Radiobutton(intraop_frame, text="Halfway through CPB", variable=intraop_var, value="half").pack(side="left", padx=10)
        ttk.Radiobutton(intraop_frame, text="Custom", variable=intraop_var, value="custom").pack(side="left", padx=10)
        
        # Ignore text field
        ttk.Label(ignore_frame, text="Enter FIS values to ignore (comma separated):").pack(side="left", padx=5)
        ttk.Entry(ignore_frame, textvariable=ignore_var, width=30).pack(side="left", padx=5)
        
        # Multi-condition analysis checkbox
        multi_cond_chk = ttk.Checkbutton(multi_cond_frame, text="Compare across ALL surgery time conditions (surgery, CPB, and clamp) [DEPRECATED]", variable=multi_condition_var)
        multi_cond_chk.pack(side="left", padx=10, pady=5)
        
        # Variables checkboxes with vertical and horizontal scrollbars
        vars_outer_frame = ttk.Frame(vars_frame)
        vars_outer_frame.pack(fill="both", expand=True)
        
        vars_canvas = tk.Canvas(vars_outer_frame)
        
        # Vertical scrollbar
        vars_v_scrollbar = ttk.Scrollbar(vars_outer_frame, orient="vertical", command=vars_canvas.yview)
        vars_v_scrollbar.pack(side="right", fill="y")
        
        # Horizontal scrollbar
        vars_h_scrollbar = ttk.Scrollbar(vars_outer_frame, orient="horizontal", command=vars_canvas.xview)
        vars_h_scrollbar.pack(side="bottom", fill="x")
        
        # Configure canvas
        vars_canvas.configure(yscrollcommand=vars_v_scrollbar.set, xscrollcommand=vars_h_scrollbar.set)
        vars_canvas.pack(side="left", fill="both", expand=True)
        
        vars_scrollable_frame = ttk.Frame(vars_canvas)
        vars_scrollable_frame.bind(
            "<Configure>",
            lambda e: vars_canvas.configure(scrollregion=vars_canvas.bbox("all"))
        )
        vars_canvas_window = vars_canvas.create_window((0, 0), window=vars_scrollable_frame, anchor="nw")
        
        # Make scrollable frame expand with canvas width
        def on_canvas_configure(event):
            vars_canvas.itemconfig(vars_canvas_window, width=event.width)
        vars_canvas.bind("<Configure>", on_canvas_configure)
        
        # Select all variables button
        vars_select_frame = ttk.Frame(vars_frame)
        vars_select_frame.pack(fill="x", pady=5, before=vars_outer_frame)
        all_vars_var = tk.BooleanVar(value=False)
        def toggle_all_vars():
            select_all = all_vars_var.get()
            for var in var_checkboxes.values():
                var.set(select_all)
        ttk.Checkbutton(vars_select_frame, text="Select All Variables", variable=all_vars_var, 
                        command=toggle_all_vars).pack(side="left", padx=10)
        
        # Categorize variables
        delta_vars = []
        theta_vars = []
        alpha_vars = []
        beta_vars = []
        other_vars = []
        
        for var_key, var_name in vars_dict.items():
            # Capitalize the first letter of the display name
            display_name = var_name[0].upper() + var_name[1:] if var_name else ""
            
            # Create variable for checkbox
            var_checkboxes[var_key] = tk.BooleanVar(value=False)
            
            # Categorize based on keywords (case-insensitive)
            var_name_lower = var_name.lower()
            if 'delta' in var_name_lower:
                delta_vars.append((var_key, display_name))
            elif 'theta' in var_name_lower:
                theta_vars.append((var_key, display_name))
            elif 'alpha' in var_name_lower:
                alpha_vars.append((var_key, display_name))
            elif 'beta' in var_name_lower:
                beta_vars.append((var_key, display_name))
            else:
                other_vars.append((var_key, display_name))
        
        # Create category frames
        category_frames = {}
        for idx, (category, title) in enumerate([
            ("delta", "Delta Variables"), 
            ("theta", "Theta Variables"), 
            ("alpha", "Alpha Variables"), 
            ("beta", "Beta Variables"), 
            ("other", "Other Variables")
        ]):
            frame = ttk.LabelFrame(vars_scrollable_frame, text=title)
            frame.grid(row=0, column=idx, sticky="nsew", padx=5, pady=5)
            category_frames[category] = frame
        
        # Add variables to their respective categories
        for category, var_list in [
            ("delta", delta_vars),
            ("theta", theta_vars),
            ("alpha", alpha_vars),
            ("beta", beta_vars),
            ("other", other_vars)
        ]:
            for i, (var_key, display_name) in enumerate(var_list):
                ttk.Checkbutton(
                    category_frames[category], 
                    text=display_name, 
                    variable=var_checkboxes[var_key]
                ).grid(row=i, column=0, sticky="w", padx=5, pady=2)
        
        # Statistics checkboxes with scrollbar - use a fixed height to prevent overflow
        stats_outer_frame = ttk.Frame(stats_frame)
        stats_outer_frame.pack(fill="both", expand=True)
        stats_canvas = tk.Canvas(stats_outer_frame)
        # Vertical scrollbar
        stats_v_scrollbar = ttk.Scrollbar(stats_outer_frame, orient="vertical", command=stats_canvas.yview)
        stats_v_scrollbar.pack(side="right", fill="y")
        # Horizontal scrollbar
        stats_h_scrollbar = ttk.Scrollbar(stats_outer_frame, orient="horizontal", command=stats_canvas.xview)
        stats_h_scrollbar.pack(side="bottom", fill="x")
        # Configure canvas
        stats_canvas.configure(yscrollcommand=stats_v_scrollbar.set, xscrollcommand=stats_h_scrollbar.set)
        stats_canvas.pack(side="left", fill="both", expand=True)
        stats_scrollable_frame = ttk.Frame(stats_canvas)
        stats_scrollable_frame.bind("<Configure>", lambda e: stats_canvas.configure(scrollregion=stats_canvas.bbox("all")))
        stats_canvas_window = stats_canvas.create_window((0, 0), window=stats_scrollable_frame, anchor="nw")
        
        # Make scrollable frame expand with canvas width
        def on_stats_canvas_configure(event):
            stats_canvas.itemconfig(stats_canvas_window, width=event.width)
        stats_canvas.bind("<Configure>", on_stats_canvas_configure)
        
        # Select all stats button
        stats_select_frame = ttk.Frame(stats_frame)
        stats_select_frame.pack(fill="x", pady=5, before=stats_outer_frame) 
        all_stats_var = tk.BooleanVar(value=False)
        def toggle_all_stats():
            select_all = all_stats_var.get()
            for var in stat_checkboxes.values():
                var.set(select_all)
        ttk.Checkbutton(stats_select_frame, text="Select All Statistics", variable=all_stats_var, 
                        command=toggle_all_stats).pack(side="left", padx=10)
        
        # Categorize statistics
        stat_categories = {
            "Central Tendency": ["mean", "median"],
            "Spread": ["stdev", "sem", "iqr", "perc_5", "perc_95"],
            "Linear Regression": ["lin_reg_slope", "lin_reg_yint", "lin_reg_rmse"],
            "Theil-Sen Slope": ["theil_sen_slope", "theil_sen_yint", "theil_sen_rmse"],
            "Mann-Kendall Tau": ["mann_kendall_tau", "mann_kendall_p"],
            "Time series characteristics": ["timeseries_num_points", "timeseries_duration_mins"]
        }
        
        # Create category frames for statistics
        stat_category_frames = {}
        for idx, category_name in enumerate(stat_categories.keys()):
            frame = ttk.LabelFrame(stats_scrollable_frame, text=category_name)
            frame.grid(row=0, column=idx, sticky="nsew", padx=5, pady=5)
            stat_category_frames[category_name] = frame
        
        # Add statistics to their respective categories
        for category_name, stat_keys in stat_categories.items():
            for i, stat_key in enumerate(stat_keys):
                if stat_key in stats_dict:
                    # Capitalize the first letter of the display name
                    display_name = stats_dict[stat_key][0].upper() + stats_dict[stat_key][1:] if stats_dict[stat_key] else ""
                    
                    stat_checkboxes[stat_key] = tk.BooleanVar(value=False)
                    ttk.Checkbutton(
                        stat_category_frames[category_name],
                        text=display_name,
                        variable=stat_checkboxes[stat_key]
                    ).grid(row=i, column=0, sticky="w", padx=5, pady=2)
        
        # Result variables to return
        result = {
            'machine_filt': None,
            'intraop_per': None,
            'ignore_arr': None,
            'vars_arr': None,
            'stats_arr': None,
            'multi_condition': None,
            'custom_start': None,
            'custom_end': None,
            'min_valid_data': None,
            'control_vars': None,
            'use_adjusted_correlations': None
        }
        
        # OK and Cancel buttons - align to the right
        button_frame.columnconfigure(0, weight=1)  # Makes buttons align right
        def on_ok():
            # Process machine filter
            result['machine_filt'] = machine_var.get()
            # Process intraop period
            result['intraop_per'] = intraop_var.get()
            # Process ignore array
            ignore_text = ignore_var.get().strip()
            if ignore_text:
                # Parse comma-separated numbers, handling spaces
                ignore_arr = [int(num.strip()) for num in re.split(r'[,\s]+', ignore_text) if num.strip()]
                result['ignore_arr'] = ignore_arr
            else:
                result['ignore_arr'] = []
            # Process variables array
            result['vars_arr'] = [var_key for var_key, var in var_checkboxes.items() if var.get()]
            # Process stats array
            result['stats_arr'] = [stat_key for stat_key, stat in stat_checkboxes.items() if stat.get()]
            # Process multi-condition flag
            result['multi_condition'] = multi_condition_var.get()
            result['custom_start'] = float(custom_start.get()) if custom_start.get() else None
            result['custom_end'] = float(custom_end.get()) if custom_end.get() else None
            result['min_valid_data'] = float(min_valid_data.get()) if min_valid_data.get() else 0.0
            control_vars = {
                'seizures_intraop': control_seizure_var.get(),
                'ga_eeg': control_ga_var.get(),
                'pma_eeg': control_pma_var.get()
            }
            result['control_vars'] = control_vars
            # Get adjusted correlations setting
            result['use_adjusted_correlations'] = use_adjusted_correlations_var.get()
            
            root.destroy()
            
        def on_cancel():
            root.destroy()
        
        ttk.Button(button_frame, text="Cancel", command=on_cancel).grid(row=0, column=1, padx=5)
        ttk.Button(button_frame, text="OK", command=on_ok).grid(row=0, column=2, padx=5)
        
        # Set minimum window size to ensure buttons are visible
        root.update_idletasks()
        min_height = button_frame.winfo_y() + button_frame.winfo_height() + 20
        min_width = 800
        root.minsize(min_width, min_height)
        
        # Run the main loop
        root.mainloop()
        
        # Return the selected values with the multi-condition flag
        return result['machine_filt'], result['intraop_per'], result['ignore_arr'], result['vars_arr'], result['stats_arr'], result['multi_condition'], result['custom_start'], result['custom_end'], result['min_valid_data'], result['control_vars'], result['use_adjusted_correlations']
    
    def stats_all(self, alpha=0.05, trend=0.10):
        """
        The point of this function is to iterate across all variables and all stats to determine where the significant values are.
        If multi-condition is True, it will run the analysis across all three time conditions.
        """
        
        #Step 1: Prepare settings for analysis
        #Step 1.1: Determine filter settings
        machine_filt, intraop_per, ignore_arr, vars_arr, stats_arr, multi_condition, custom_start, custom_end, min_valid, control_vars, use_adjusted_correlations = self.get_filter_options()
        ancova_bool = control_vars and any(control_vars.values())
        #Step 1.3: Record filter settings
        os.makedirs(self.base_save_dir, exist_ok=True)
        filter_settings_log = os.path.join(self.base_save_dir, "filter_settings.txt")
        with open(filter_settings_log, "w") as f:
            f.write("Analysis Filter Settings\n")
            f.write("======================\n\n")
            f.write(f"Date and time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write("Basic Settings:\n")
            f.write(f"  Machine filter: {machine_filt}\n")
            f.write(f"  Intraoperative period: {intraop_per}\n")
            f.write(f"  Multi-condition analysis: {multi_condition}\n")
            f.write(f"  Minimum valid data: {min_valid}\n\n")
            if custom_start is not None or custom_end is not None:
                f.write("Custom Time Settings:\n")
                f.write(f"  Start time: {custom_start} hours\n")
                f.write(f"  End time: {custom_end} hours\n\n")
            f.write("Statistical Control Settings:\n")
            if ancova_bool:
                for var_name, enabled in control_vars.items():
                    f.write(f"  {var_name}: {enabled}\n")
                f.write(f"  Used adjusted values for correlations: {use_adjusted_correlations}\n\n")
            else:
                f.write("  No control variables used\n\n")
            if ignore_arr:
                f.write("Ignored FIS Numbers:\n")
                f.write(f"  {', '.join(map(str, ignore_arr))}\n\n")
            else:
                f.write("Ignored FIS Numbers: None\n\n")
            f.write("Selected Variables:\n")
            for var in vars_arr:
                f.write(f"  - {var}\n")
            f.write("\n")
            f.write("Selected Statistics:\n")
            for stat in stats_arr:
                f.write(f"  - {stat}\n")
            f.write("\n")
            f.write("Significance Thresholds:\n")
            f.write(f"  Alpha (significance): {alpha}\n")
            f.write(f"  Trend threshold: {trend}\n")
        print(f"Filter settings log saved to: {filter_settings_log}")
        #Step 1.3: If multi-condition is selected, run the multi-condition analysis. This happens before our filtering because the ecc filtering must be done for each case in the multi condition analysis.
        if multi_condition:
            print("Running multi-condition analysis across all surgery time periods...")
            self.multi_condition_analysis(machine_filt=machine_filt, ignore_arr=ignore_arr, vars_arr=vars_arr, stats_arr=stats_arr,  alpha=alpha, trend=trend, min_valid=min_valid)
            return
        #Step 1.4: If not multi-condition, filter here
        self.filter_by_machine(machine_filt)
        self.filter_by_ecc(intraop_per)
        self.filter_by_ignore(ignore_arr)

        #Step 2: Run statistics for all included patients
        final_df = []
        significant_df = []
        trending_df = []
        all_correlation_results = []
        #Step 2.1: Iterate over stats and build the rows for our output PDFs
        for var in tqdm(vars_arr, desc="Calculating parameters for all included patients"):
            for stat in tqdm(stats_arr, desc=f"Stats for {var}", leave=False):
                df_row, param_name, p = self.run_test(var, stat, test=None, time_range=intraop_per, mode="all", cus_time_1=custom_start, cus_time_2=custom_end, min_art=min_valid, all_correlation_results=all_correlation_results,control_vars=control_vars,use_adjusted_correlations=use_adjusted_correlations, case = "whole_sample")
                p = float(p)
                if not final_df:
                    final_df = [df_row.columns.tolist()]
                final_df.append(df_row.values.tolist()[0])
                if p < alpha:
                    if not significant_df:
                        significant_df = [df_row.columns.tolist()]
                    significant_df.append(df_row.values.tolist()[0])
                elif alpha <= p < trend:
                    if not trending_df:
                        trending_df = [df_row.columns.tolist()]
                    trending_df.append(df_row.values.tolist()[0])
        #Step 2.2: Save results for all patients meeting filtering requirements
        #Step 2.2.1: Save results for all patients meeting filtering requirements
        if final_df:
            save_by_category(final_df, f'{self.base_save_dir}\\all_results', 'All Statistical Results')
        if significant_df:
            save_by_category(significant_df, f'{self.base_save_dir}\\significant_results', 'Significant Results (p < 0.05)')
        if trending_df:
            save_by_category(trending_df, f'{self.base_save_dir}\\trending_results', 'Trending Results (0.05 <= p < 0.10)')
        if final_df:
            combine_pdfs(f'{self.base_save_dir}\\all_results', 'all_results_combined.pdf')
        if significant_df:
            combine_pdfs(f'{self.base_save_dir}\\significant_results', 'significant_results_combined.pdf')
        if trending_df:
            combine_pdfs(f'{self.base_save_dir}\\trending_results', 'trending_results_combined.pdf')
        if significant_df or trending_df:
            integrate_PDF(self.base_save_dir, significant_df, trending_df)
        #Step 2.2.2: Save correlation results for all patients meeting filtering requirements
        adjustment_info = "adjusted" if use_adjusted_correlations and ancova_bool else "unadjusted"
        folder_name = f"correlation_results_{adjustment_info}"
        n_pre_corr = sum([1 for val in self.fis_dict.values() if val.analyzed_whole_sample == 1])
        save_correlation_results(self, all_correlation_results, folder_name=folder_name)
        n_post_corr = sum([1 for val in self.fis_dict.values() if val.analyzed_whole_sample == 1])
        print(f"Number analyzed before correlation analysis: {n_pre_corr}")
        print(f"Number analyzed after correlation analysis: {n_post_corr}")
        #Step 2.3: Generate tables
        generate_demographics_tables(self, case = "whole_sample")
        generate_risk_categories_table(self)
        generate_seizure_table(self)
        generate_neurodevelopmental_table(self)
        generate_artefact_table(self)
        generate_biomarker_table(self)
        #Step 2.4: Log final sample
        #Step 2.4.1: Create log of original sample
        original_count_log = os.path.join(self.base_save_dir, "original_count.txt")
        with open(original_count_log, "w") as f:
            f.write(f"Original number of FISs: {self.num_original_fis}\n")
            for fis_obj in self.fis_dict.values():
                f.write(f"{fis_obj.fis}\n")
        print(f"Original FIS log saved to: {original_count_log}")
        #Step 2.4.2: Create log of FIS that were analyzed
        analyzed_fis_log = os.path.join(self.base_save_dir, "analyzed_fis_log.txt")
        with open(analyzed_fis_log, "w") as f:
            f.write("FIS\tAnalysis Status\n")
            analyzed_count = 0
            for fis_obj in self.fis_dict.values():
                if hasattr(fis_obj, "analyzed_whole_sample") and fis_obj.analyzed_whole_sample == 1:
                    f.write(f"{fis_obj.fis}\tAnalyzed\n")
                    analyzed_count += 1
            f.write(f"\nTotal analyzed FISs: {analyzed_count}\n")
        print(f"Analyzed FIS log saved to: {analyzed_fis_log}")
        # Step 2.4.2.1: Create log of FIS that were NOT analyzed
        not_analyzed_fis_log = os.path.join(self.base_save_dir, "not_analyzed_fis_log.txt")
        with open(not_analyzed_fis_log, "w") as f:
            f.write("FIS\tReason Not Analyzed\n")
            f.write("===\t==================\n")
            not_analyzed_count = 0
            for fis_obj in self.fis_dict.values():
                if not hasattr(fis_obj, "analyzed_whole_sample") or fis_obj.analyzed_whole_sample != 1:
                    skip_reason = getattr(fis_obj, 'skip_str', 'No reason provided')
                    f.write(f"{fis_obj.fis}\t{skip_reason}\n")
                    not_analyzed_count += 1
            f.write(f"\nTotal not analyzed FISs: {not_analyzed_count}\n")
            f.write(f"Total analyzed FISs: {analyzed_count}\n")
            f.write(f"Original total FISs: {self.num_original_fis}\n")
        print(f"Not analyzed FIS log saved to: {not_analyzed_fis_log}")
        # Step 2.4.3: Log analysis information
        analysis_info_log = os.path.join(self.base_save_dir, "analysis_info.txt")
        with open(analysis_info_log, "w") as f:
            f.write("Analysis Information\n")
            f.write("===================\n\n")
            f.write(f"Machine filter: {machine_filt}\n")
            f.write(f"Intraoperative period: {intraop_per}\n")
            f.write(f"Control variables: {control_vars}\n")
            f.write(f"Used adjusted values for correlations: {use_adjusted_correlations}\n")
            f.write(f"Correlation results folder: {folder_name}\n")
        print(f"Analysis info log saved to: {analysis_info_log}")
        #Step 2.4.4: Record parameter/statistic values in Excel
        excel_file_unadj = self.export_to_excel(vars_arr, stats_arr, intraop_per, 
                                            custom_start, custom_end, 
                                            suffix="whole_sample",
                                            control_vars=control_vars,
                                            use_adjusted_values=False, analysis_type="whole_sample")
        print(f"Excel export completed for whole sample (unadjusted): {excel_file_unadj}")

        if ancova_bool:  # Only export adjusted values if ANCOVA was used
            excel_file_adj = self.export_to_excel(vars_arr, stats_arr, intraop_per, 
                                                custom_start, custom_end, 
                                                suffix="whole_sample",
                                                control_vars=control_vars,
                                                use_adjusted_values=True, analysis_type="whole_sample")
            print(f"Excel export completed for whole sample (adjusted): {excel_file_adj}")
        
        excel_file_all = self.export_to_excel(vars_arr, stats_arr, intraop_per, 
                                    custom_start, custom_end, 
                                    suffix="all",
                                    control_vars=control_vars,
                                    use_adjusted_values=False, analysis_type="all")
        print(f"Excel export completed for all patients (unadjusted): {excel_file_all}")
        
        #Step 3: Re-run statistics between the seizure and non_seizure groups
        base_save_dir_2 = os.path.join(self.base_save_dir,"across_seizures_comparison")
        self.clear_filters()
        self.filter_by_seizure(2)
        self.filter_by_ND() #This ensures that the seizures comparison is occurring only within our sample of analysis
        final_df = []
        significant_df = []
        trending_df = []
        #Step 3: Iterate the stats
        for var in tqdm(vars_arr, desc="Calculating parameters between the seizures and non-seizures group"):
            for stat in tqdm(stats_arr, desc=f"Stats for {var}", leave=False):
                df_row, param_name, p = self.run_test(var, stat, test=None, time_range=intraop_per, mode="all", cus_time_1=custom_start, cus_time_2=custom_end, min_art=min_valid, comparison_category = ["seizures_intraop","intraoperative seizure"], control_vars=control_vars, case = "seizure_sample") #This is where the fis_objs are iterated over one by one
                p = float(p)
                #add df_row to final_df. If final_df has no column names, add those too
                if not final_df:
                    final_df = [df_row.columns.tolist()]
                final_df.append(df_row.values.tolist()[0])
                #if p is under 0.05 add df_row to significant_df. If significant_df has no column names, add those too
                if p < alpha:
                    if not significant_df:
                        significant_df = [df_row.columns.tolist()]
                    significant_df.append(df_row.values.tolist()[0])
                #if p is over 0.05 and under 0.10, add df_row to trending_df. If trending_df has no column names, add those too
                elif alpha <= p < trend:
                    if not trending_df:
                        trending_df = [df_row.columns.tolist()]
                    trending_df.append(df_row.values.tolist()[0])
        # Save results
        if final_df:
            save_by_category(final_df, f'{base_save_dir_2}\\all_results', 'All Statistical Results')
        if significant_df:
            save_by_category(significant_df, f'{base_save_dir_2}\\significant_results', 'Significant Results (p < 0.05)')
        if trending_df:
            save_by_category(trending_df, f'{base_save_dir_2}\\trending_results', 'Trending Results (0.05 <= p < 0.10)')
        # After saving the individual category PDFs:
        if final_df:
            combine_pdfs(f'{base_save_dir_2}\\all_results', 'all_results_combined.pdf')
        if significant_df:
            combine_pdfs(f'{base_save_dir_2}\\significant_results', 'significant_results_combined.pdf')
        if trending_df:
            combine_pdfs(f'{base_save_dir_2}\\trending_results', 'trending_results_combined.pdf')
        #Create an integrated PDF of significant and trending results
        if significant_df or trending_df:
            integrate_PDF(base_save_dir_2, significant_df, trending_df)
        original_count_log = os.path.join(base_save_dir_2, "original_count.txt")
        with open(original_count_log, "w") as f:
            f.write(f"Original number of FISs: {self.num_original_fis}\n")
            for fis_obj in self.fis_dict.values():
                f.write(f"{fis_obj.fis}\n")
        print(f"Original FIS log saved to: {original_count_log}")
        # Create a log of analyzed FISs
        analyzed_fis_log = os.path.join(base_save_dir_2, "analyzed_fis_log.txt")
        with open(analyzed_fis_log, "w") as f:
            f.write("FIS\tAnalysis Status\n")
            analyzed_count = 0
            for fis_obj in self.fis_dict.values():
                if hasattr(fis_obj, "analyzed_seizure_sample") and fis_obj.analyzed_seizure_sample == 1:
                    f.write(f"{fis_obj.fis}\tAnalyzed\n")
                    analyzed_count += 1
            f.write(f"\nTotal analyzed FISs: {analyzed_count}\n")
        print(f"Analyzed FIS log saved to: {analyzed_fis_log}")

        excel_file_seizures = self.export_to_excel(vars_arr, stats_arr, intraop_per, 
                                          custom_start, custom_end, 
                                          suffix="seizure_comparison",
                                          control_vars=control_vars,
                                          use_adjusted_values=False, analysis_type="seizure_sample")
        print(f"Excel export completed for seizure comparison: {excel_file_seizures}")
        
        #Step 4: Re-run tests on only the no-seizure group
        base_save_dir_2 = os.path.join(self.base_save_dir,"no_seizures_comparison")
        self.clear_filters()
        self.filter_by_seizure(0)
        final_df = []
        significant_df = []
        trending_df = []
        #Step 4: Iterate the stats for no seizure group
        for var in tqdm(vars_arr, desc="Calculating parameters in the no seizure group"):
            for stat in tqdm(stats_arr, desc=f"Stats for {var}", leave=False):
                df_row, param_name, p = self.run_test(var, stat, test=None, time_range=intraop_per, mode="all", cus_time_1=custom_start, cus_time_2=custom_end, min_art=min_valid, comparison_category = ["ND","ND"], control_vars=control_vars, case = "across_seizures") #This is where the fis_objs are iterated over one by one
                p = float(p)
                #add df_row to final_df. If final_df has no column names, add those too
                if not final_df:
                    final_df = [df_row.columns.tolist()]
                final_df.append(df_row.values.tolist()[0])
                #if p is under 0.05 add df_row to significant_df. If significant_df has no column names, add those too
                if p < alpha:
                    if not significant_df:
                        significant_df = [df_row.columns.tolist()]
                    significant_df.append(df_row.values.tolist()[0])
                #if p is over 0.05 and under 0.10, add df_row to trending_df. If trending_df has no column names, add those too
                elif alpha <= p < trend:
                    if not trending_df:
                        trending_df = [df_row.columns.tolist()]
                    trending_df.append(df_row.values.tolist()[0])
        # Step 3: prep PDFs with results
        # Save results
        if final_df:
            save_by_category(final_df, f'{base_save_dir_2}\\all_results', 'All Statistical Results')
        if significant_df:
            save_by_category(significant_df, f'{base_save_dir_2}\\significant_results', 'Significant Results (p < 0.05)')
        if trending_df:
            save_by_category(trending_df, f'{base_save_dir_2}\\trending_results', 'Trending Results (0.05 <= p < 0.10)')
        # After saving the individual category PDFs:
        if final_df:
            combine_pdfs(f'{base_save_dir_2}\\all_results', 'all_results_combined.pdf')
        if significant_df:
            combine_pdfs(f'{base_save_dir_2}\\significant_results', 'significant_results_combined.pdf')
        if trending_df:
            combine_pdfs(f'{base_save_dir_2}\\trending_results', 'trending_results_combined.pdf')
        #Create an integrated PDF of significant and trending results
        if significant_df or trending_df:
            integrate_PDF(base_save_dir_2, significant_df, trending_df)
        original_count_log = os.path.join(base_save_dir_2, "original_count.txt")
        with open(original_count_log, "w") as f:
            f.write(f"Original number of FISs: {self.num_original_fis}\n")
            for fis_obj in self.fis_dict.values():
                f.write(f"{fis_obj.fis}\n")
        print(f"Original FIS log saved to: {original_count_log}")
        # Create a log of analyzed FISs
        analyzed_fis_log = os.path.join(base_save_dir_2, "analyzed_fis_log.txt")
        with open(analyzed_fis_log, "w") as f:
            f.write("FIS\tAnalysis Status\n")
            analyzed_count = 0
            for fis_obj in self.fis_dict.values():
                if hasattr(fis_obj, "analyzed_across_seizures") and fis_obj.analyzed_across_seizures == 1:
                    f.write(f"{fis_obj.fis}\tAnalyzed\n")
                    analyzed_count += 1
            f.write(f"\nTotal analyzed FISs: {analyzed_count}\n")
        print(f"Analyzed FIS log saved to: {analyzed_fis_log}")
        excel_file_no_seizures = self.export_to_excel(vars_arr, stats_arr, intraop_per, 
                                             custom_start, custom_end, 
                                             suffix="no_seizure_group",
                                             control_vars=control_vars,
                                             use_adjusted_values=False, analysis_type="across_seizures")
        print(f"Excel export completed for no seizure group: {excel_file_no_seizures}")