r"""
This code is intended to be used after both OBM_stitching_cutting.py and Nicolet_stitching_cutting.py have been run (as it relies on some of the output Excel sheets from those functions). The purpose of this file is to combine the signals from multiple alex_raw_efficient objects that are in the same temporal window, such that these signals can be analyzed all as one signal.

NB: in a given window of interest (e.g., 0 to 24 hrs postop), we may have multiple files in "D:\EEG_Pipeline_Continuous\0_raw_edf_files\downloaded_nicolet_files" (or .../OBM_files) that contain data in that window. Therefore, the big purpose of this code is to find separate alex_raw_efficient objects whose data is in the same chunk of interest, and then add in the data from this chunk so it is a continuous data stream (with zeros added as needed) for each window. Therefore, the output will be one alex_raw_efficient object (saved as a .pkl file to D:\EEG_Pipeline_Continuous\1_stitched_cut\final_stitched_cut_efficient) with the chunk it belongs to specified in the file name.
"""
import pandas as pd
import os
import sys
import shutil
current_dir = os.path.dirname(os.path.abspath(__file__))
utils_dir = os.path.join(current_dir, '..', 'utils')
sys.path.append(utils_dir)
from file_importing import import_pkl
import ast
import numpy as np
import copy
import re
import ast
from datetime import datetime
import tqdm as tqdm

def record_error(error_string):
    """
    Helper function that builds an error log for this fxn.

    Inputs:
    - error_string (str): the string we want added to the error log
    """
    # Define the log file path
    log_file_path = "error_log_file_combiner.txt"

    # Write the error message to the log file
    with open(log_file_path, "a") as log_file:
        # Write log initialization with today's date
        today_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # Format date and time
        log_file.write(f"Date of error: {today_date}\n")
        
        # Write the error string to the file
        log_file.write(f"- {error_string}\n\n")

def blend(save_base_dir, suffix, windows, clear = True):
    """
    The purpose of this code is to stitch together multiple Nicolet alex_raw_efficient.pkl files in the same chunk into one that contains all the data (and only the data) from that chunk.

    Inputs:
    - read_file (str): path to Excel file containing the FIS and chunk data
    - suffix (str): suffix that we add when saving the file (to help us know nicolet ["N"] or OBM ["OBM"])
    """
    output_path = os.path.join(save_base_dir, "final_stitched_cut_efficient")
    if suffix == "N":
        data_sheet_path = os.path.join(save_base_dir, "nicolet_raw_file_timestamps.xlsx")
    elif suffix == "O":
        data_sheet_path = os.path.join(save_base_dir, "OBM_raw_file_timestamps.xlsx")
    df = pd.read_excel(data_sheet_path)
    df["fis"] = df["fis"].apply(lambda x: int(x) if isinstance(x, str) and x.isdigit() else x)

    #Step 0: Clear if necessary
    os.makedirs(output_path, exist_ok=True)
    if clear:
        for item in os.listdir(output_path):
            item_path = os.path.join(output_path, item)
            if os.path.isfile(item_path):
                os.remove(item_path)
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
    
    #Step 1: For each combination of FIS and chunk, find all files of that given FIS in the time chunk of interest.
    for fis in df["fis"].unique():
        #Find all files of a given FIS
        fis_subset = df[(df["fis"] == fis) | (df["fis"].astype(str) == str(fis))]
        #For each row in fis_subset, examine the values in the columns start_time and end_time. If they fall into the time range denoted by

        #For that FIS, find all files in a given chunk
        for chunk in fis_subset["chunk"].unique():
            if chunk == "-": #Hard-coded error case
                continue
            chunk_subset = fis_subset[fis_subset["chunk"] == chunk]
            #Get the paths of all these files
            stitched_paths = [entry for entry in chunk_subset["stitched_path"]]
            chunk = list(map(float, ast.literal_eval(chunk)))
            chunk_list = [val for val in chunk]
            chunk = [val * 3600 for val in chunk]
            base_data = []
            base_times_from_zero = []
            base_times_postop = []
            srates = []

            #Step 2: For each FIS-chunk combo of interest, get the data and time arrays
            print(f"\nNow examining FIS{fis} at chunk {[round(val/3600,0) for val in chunk]}")
            for i, path in enumerate(stitched_paths):
                if suffix == "O":
                    file_base = os.path.split(path)[0]
                    file_name = os.path.split(path)[1]
                    file_name = file_name.split("_")[0]
                    file_name = f"{file_name}_stitched_cut.pkl"
                    path = os.path.join(file_base,file_name)
                raw_obj = import_pkl(full_file_path=path, parts_or_path="path")
                base_data.append(raw_obj.EEG_data)
                base_times_from_zero.append(raw_obj.EEG_times_from_zero)
                base_times_postop.append(raw_obj.EEG_times_from_analyze_end_time)
                srates.append(raw_obj.sampling_rate)
                if i == len(stitched_paths) - 1:  # Check if this is the last iteration
                    fin_raw = copy.deepcopy(raw_obj)
                del raw_obj

            fis_pattern = re.search(r'(?i)FIS\d+', fin_raw.file_name)
            fis_id = fis_pattern.group(0) if fis_pattern else fin_raw.file_name.split('_')[0]
            fin_path = output_path + "\\" + fis_id + f"_{int(chunk_list[0])}_{int(chunk_list[1])}_{suffix}.pkl"


            #Step 3: Order base_times_postop and apply the same ordering to base_data and base_times_from_zero
            sorted_indices = sorted(range(len(base_times_postop)), key=lambda i: base_times_postop[i][0]) #sorted_indices is an array where each entry is the index of the entries in base_times_postop in ascending order ([0,2,1] means that item 0 of base_times_postop occurred first, then item 2, then item 1)
            base_times_postop = [base_times_postop[i] for i in sorted_indices] #This step simply reorders base_times_postop to be in ascending order (by first entry)
            base_times_from_zero = [base_times_from_zero[i] for i in sorted_indices]
            base_data = [base_data[i] for i in sorted_indices]
            #analyze_end_time = base_times_postop[0][0] #analyze_end_time is our first postop time point
            # View the postop time ranges easily
            postop_time_ranges = {(arr[0], arr[-1]) for arr in base_times_postop}
            print(f"Postop time ranges in the original signal: {postop_time_ranges}")

            #Check for and handle overlaps across files
            def overlap_check():
                nonlocal base_times_postop, base_times_from_zero, base_data

                overlap_found = True
                while overlap_found:
                    i = 0
                    overlap_found = False
                    postop_time_ranges = [(arr[0], arr[-1]) for arr in base_times_postop]
                    print(len(postop_time_ranges))
                    #This loop checks and fixes all, but we should break
                    for i in range(len(postop_time_ranges)-1):
                        curr_arr_last_time = postop_time_ranges[i][1]
                        next_arr_first_time = postop_time_ranges[i+1][0]
                        if curr_arr_last_time > next_arr_first_time:  # In this case, we have an overlap that we need to fix
                            # Find the indices where the values in base_times_postop[i+1] are less than the last value of base_times_postop[i]
                            overlap_indices = [j for j, time in enumerate(base_times_postop[i+1]) if time <= curr_arr_last_time]
                            if not len(base_data[i+1]) == len(base_times_postop[i+1]) and len(base_data[i+1]) == len(base_times_from_zero[i+1]):
                                raise ValueError("The data individual data and time segments for a given file should have the same length, but they do not.")
                            # Eliminate all items at these indices in base_data[i+1], base_times_postop[i+1] and base_times_from_zero[i+1]
                            overlap_set = set(overlap_indices)
                            # Create a filtered copy of each list in one pass
                            keep_indices = [j for j in range(len(base_data[i+1])) if j not in overlap_set]
                            base_data[i+1] = [base_data[i+1][j] for j in keep_indices]
                            base_times_postop[i+1] = [base_times_postop[i+1][j] for j in keep_indices]
                            base_times_from_zero[i+1] = [base_times_from_zero[i+1][j] for j in keep_indices]

                            base_times_postop = [arr for arr in base_times_postop if not len(arr)==0] #Handles case where one recording is totally inside another
                            postop_time_ranges = [(arr[0], arr[-1]) for arr in base_times_postop]
                            curr_arr_last_time = postop_time_ranges[i][1]
                            next_arr_first_time = postop_time_ranges[i+1][0]
                            print(f"{curr_arr_last_time} sec should be before {next_arr_first_time} sec.")
                            if not curr_arr_last_time < next_arr_first_time:
                                raise ValueError("The current last time should be before the subsequent first time.") 
                            overlap_found = True
                            break #Takes us back to while loop to re-run with the updated variables

            overlap_check()
            
            #Step 4: Parse down values to be only within the times of interest
            filtered_data = []
            filtered_times_from_zero = []
            filtered_times_postop = []
            for data, times_zero, times_postop in zip(base_data, base_times_from_zero, base_times_postop):
                # Find indices where times_postop are within the chunk
                indices = [i for i, t in enumerate(times_postop) if chunk[0] <= t <= chunk[1]]
                if indices:
                    filtered_data.append(data[indices[0]:indices[-1] + 1])
                    filtered_times_from_zero.append(times_zero[indices[0]:indices[-1] + 1])
                    filtered_times_postop.append(times_postop[indices[0]:indices[-1] + 1])
                else:
                    error_str = f"Error in the file {fin_path}: No data found in the chunk of interest despite the chunk being mentioned in nicolet_fis_and_chunk.xlsx."
                    record_error(error_str)
                    continue
            #View filtered postop time ranges
            postop_time_ranges = [(arr[0], arr[-1]) for arr in filtered_times_postop]
            print(f"Postop time ranges after filtering by chunk: {postop_time_ranges}")
            if any(postop_time_ranges[i][1] > postop_time_ranges[i + 1][0] for i in range(len(postop_time_ranges) - 1)):
                error_str = f"Error in the file {fin_path}: Postop time ranges ({postop_time_ranges}) are not in ascending order."
                record_error(error_str)
                continue
            
            if len(set(srates)) > 1:
                error_str = f"Error in the file {fin_path}: Inconsistent sampling rates detected: {srates}."
                record_error(error_str)
                continue
            #Step 5: Combine all arrays in filtered_data into one continuous array, adding zeros as needed in between arrays
            combined_data = []
            srate = srates[0]
            for i in range(len(filtered_data) - 1):
                combined_data.extend(filtered_data[i])
                time_gap = filtered_times_postop[i + 1][0] - filtered_times_postop[i][-1] #time_gap is in seconds
                gap_samples = int(time_gap * srate)
                if gap_samples > 0:
                    combined_data.extend([0] * gap_samples)
            combined_data.extend(filtered_data[-1])
            # Now that combined_data is built, build times_from_zero array starting at zero and incrementing by 1/srate
            times_from_zero = np.linspace(0, (len(combined_data) - 1) / srate, len(combined_data))
            # Build times_postop array starting at analyze_end_time and incrementing by 1/srate sec
            times_postop = np.linspace(postop_time_ranges[0][0], filtered_times_postop[-1][-1], num=len(combined_data))
            # Convert combined_data to a numpy array
            combined_data = np.array(combined_data)

            if not (len(times_from_zero) == len(times_postop) == len(combined_data)):
                error_str = f"Error in the file {fin_path}: Length mismatch: times_from_zero ({len(times_from_zero)}), times_postop ({len(times_postop)}), combined_data ({len(combined_data)})"
                record_error(error_str)
                continue

            #Step 6: Update object metadata and resave file
            fin_raw.windows = [chunk]
            fin_raw.EEG_data = combined_data
            fin_raw.EEG_times_from_analyze_end_time = times_postop
            fin_raw.EEG_times_from_zero = times_from_zero
            chunk_list = [int(val/3600) for val in chunk]  # Convert string to list
            new_path = re.sub(r'(.*\\)FIS(\d+)[^\d].*\.pkl', r'\1FIS\2.pkl', fin_raw.file_path_current)
            fin_raw.file_path_current = new_path
            base_path = output_path + "\\" + fin_raw.file_name + f"_{chunk_list[0]}_{chunk_list[1]}_"
            matching_files = [f for f in os.listdir(output_path) if f.startswith(os.path.basename(base_path)) and f.endswith('.pkl')]
            #Checks if this FIS/path combo has already been saved (which can occur if we have both an OBM and a Nicolet recording). In these cases, the file with more data is retained
            if matching_files:
                temp_raw = import_pkl(full_file_path=fin_path, parts_or_path="path")
                # Count non-zero values in both EEG_data arrays
                temp_nonzero = np.count_nonzero(temp_raw.EEG_data)
                fin_nonzero = np.count_nonzero(fin_raw.EEG_data)
                if temp_nonzero >= fin_nonzero:
                    pass
                else:
                    fin_raw.save_as_easy(output_path, f"_{chunk_list[0]}_{chunk_list[1]}_{suffix}")
            else:
                fin_raw.save_as_easy(output_path, f"_{chunk_list[0]}_{chunk_list[1]}_{suffix}")


def combine_files(save_base_dir, windows):
    blend(save_base_dir, "N", windows, clear = True)
    try:
        blend(save_base_dir, "O", windows, clear = False)
    except FileNotFoundError:
        print("No OBM files seem to have been analyzed.")