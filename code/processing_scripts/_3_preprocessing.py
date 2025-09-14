"""
EEG Signal Preprocessing and Frequency Band Filtering Pipeline

This module applies signal preprocessing to artifact-rejected EEG data,
generating multiple frequency-specific representations for spectral analysis. The pipeline
implements clean-run segmentation to ensure filtering is only applied to continuous
artifact-free data segments.

Purpose:
--------
Transforms artifact-rejected alex_raw_efficient objects into multiple filtered versions
representing different EEG frequency bands.

Core Processing Pipeline:
------------------------
1. Clean Run Detection: Identify continuous segments of artifact-free data
2. Minimum Segment Filtering: Convert segments shorter than threshold to NaN
3. Multi-Band Filtering: Apply frequency-specific Butterworth filters
4. Notch Filtering: Remove power line interference (50 or 60 Hz, depending on continent)
5. Padding Strategy: Use reflection padding to minimize edge artifacts

Clean Run Segmentation:
----------------------
- Identifies continuous runs of valid (non-NaN) data
- Rejects (i.e., marks with NaN) segments shorter than min_seg_length to avoid filter artifacts
- Preserves temporal structure by maintaining NaN regions
- Creates clean_runs and clean_runs_lengths attributes for tracking

Filtering Implementation:
------------------------
- Butterworth filters: Zero-phase filtering using filtfilt()
- Reflection padding: 3 cycles of lowest frequency (from Mike X Cohen, Analyzing Neural Time Series Data)
- Quality control: Length validation to ensure no data loss
- Segmented processing: Filters applied only to clean runs

Output Structure:
----------------
Each input object generates 6 filtered versions:
1. Notch filtered (full bandpass spectrum, no power line interference)
2. Delta band filtered  
3. Theta band filtered
4. Alpha band filtered
5. Beta band filtered
6. Original (trimmed but unfiltered)

Configuration Dependencies:
--------------------------
Required variables from config.get_vars_dict():
- min_seg_length: Minimum segment duration (seconds) for filtering
- bandpass_order: Butterworth filter order
- notch_frequency: Power line frequency (50 or 60 Hz)
- notch_q: Quality factor for notch filter
- full_range: [low, high] frequencies for broadband
- delta_range: [low, high] frequencies for delta band
- theta_range: [low, high] frequencies for theta band  
- alpha_range: [low, high] frequencies for alpha band
- beta_range: [low, high] frequencies for beta band

File Organization:
-----------------
Input: artefact_rejected_files/ directory
Output directories:
- preprocessed_files/: Notch/bandpass filtered objects
- preprocessed_files_delta/: Delta band objects
- preprocessed_files_theta/: Theta band objects  
- preprocessed_files_alpha/: Alpha band objects
- preprocessed_files_beta/: Beta band objects
- preprocessed_files_original/: Trimmed originals

Error Handling:
--------------
- Validates filter output lengths match input lengths
- Checks that notch filtering actually modifies data
- Handles edge cases in run detection
- Graceful handling of completely artifact files

Performance Considerations:
--------------------------
- Memory management: Explicit deletion of processed objects
- Segmented filtering: Only processes valid data regions
- Deepcopy usage: Prevents cross-contamination between frequency bands
- Reflection padding: Minimizes computational overhead vs. other padding methods
"""

import numpy as np
from scipy.signal import butter, iirnotch, filtfilt
import os
import sys
current_dir = os.path.dirname(os.path.abspath(__file__))
utils_dir = os.path.join(current_dir, '..', 'utils')
sys.path.append(utils_dir)
from utils import config
vars_dict = config.get_vars_dict()
path = r'D:\utils'
sys.path.append(path)
from file_importing import import_pkl
import shutil
from copy import deepcopy

def bandpass_filter_signal(my_raw, passband):
    """
    Apply Butterworth bandpass filter to EEG signal with artifact-aware processing.
    
    Implements filtering that respects artifact boundaries and uses
    reflection padding to minimize edge effects. Only processes continuous clean
    runs identified in the clean_runs attribute.
    
    Parameters:
    -----------
    my_raw : alex_raw object
        EEG object containing signal data and clean run information
        Must have attributes: EEG_data, sampling_rate, clean_runs, clean_runs_lengths
    passband : list or array [low, high]
        Frequency band boundaries in Hz
        Example: [0.5, 4] for delta band, [8, 13] for alpha band
        
    Returns:
    --------
    my_raw : alex_raw object
        Input object with EEG_data modified in-place
        Clean runs filtered, NaN regions preserved unchanged
        
    Algorithm Details:
    -----------------
    1. Butterworth Filter Design:
       - Zero-phase filtering using filtfilt() to prevent phase distortion
       - Order specified by bandpass_order configuration parameter
       - Normalized frequencies relative to Nyquist frequency
       
    2. Padding Strategy (Mike X Cohen standard):
       - Reflection padding: 3 cycles of lowest frequency
       - Formula: pad_time = 3 / low_frequency seconds
       - Minimizes edge artifacts while preserving signal characteristics
       
    3. Segmented Processing:
       - Processes only clean_runs segments (continuous artifact-free data)
       - Preserves NaN regions to maintain temporal structure
       - Independent filtering of each clean segment prevents cross-contamination
       
    4. Quality Control:
       - Validates output length matches input length
       - Ensures no data loss during filtering process
    """
    sampling_rate = my_raw.sampling_rate
    lowpass = passband[0]
    highpass = passband[1]
    order = vars_dict["bandpass_order"]

    def butter_bandpass(lowcut=lowpass, highcut=highpass, fs=sampling_rate, order=order):
        nyquist = 0.5 * fs  # Nyquist frequency
        low = lowcut / nyquist
        high = highcut / nyquist
        b, a = butter(order, [low, high], btype='band')
        return b, a

    def butter_bandpass_filter(data=my_raw.EEG_data, lowcut=lowpass, highcut=highpass, fs=sampling_rate, order=order):
        b, a = butter_bandpass(lowcut, highcut, fs, order=order)
        time_to_pad = 1/lowcut*3 #Per Analyzing Neural Time Series Data by Mike X Cohen, use 3 cycles of the lowest frequency to pad (page 77); this is time_to_pad in SECONDS
        n_to_pad = int(time_to_pad*fs) #Convert time into number of values
        data_padded = np.pad(data, n_to_pad, mode='reflect')
        y = filtfilt(b, a, data_padded) #Filtfilt to avoid phase disruptions
        y = y[n_to_pad:-n_to_pad] #Remove padding
        return y
    
    time_to_pad = 1/passband[0]*3 #Per Analyzing Neural Time Series Data by Mike X Cohen, use 3 cycles of the lowest frequency to pad (page 77); this is time_to_pad in SECONDS
    n_to_pad = int(time_to_pad*my_raw.sampling_rate) #Convert time into number of values
    
    for start, length in zip(my_raw.clean_runs, my_raw.clean_runs_lengths):
        data_to_filter = my_raw.EEG_data[start:start+length]
        data_padded = np.pad(data_to_filter, n_to_pad, mode='reflect')
        filtered_segment = butter_bandpass_filter(data=data_padded)
        filtered_segment = filtered_segment[n_to_pad:-n_to_pad] #Remove padding
        if len(data_to_filter) != len(filtered_segment):
            raise ValueError("The length of the EEG data segment and the filtered segment should be the same.")
        my_raw.EEG_data[start:start+length] = filtered_segment
    
    my_raw.add_change_log(f"Data bandpass filtered to the range {passband} Hz with a Butterworth filter of order {order} and a sampling rate of {sampling_rate} Hz.")
    return my_raw

def notch_filter_signal(my_raw):
    """
    Code to create and apply a notch filter to the my_raw object of interest.
    
    Inputs:
    - my_raw (alex_raw object): object whose data we want to notch filter

    Outputs:
    - my_raw (alex_raw object): returns the alex_raw object with the updated EEG_data array
    """
    notch_frequency = vars_dict["notch_frequency"]
    q_factor = vars_dict["notch_q"]
    
    def iir_notch_filter(data, notch_freq, quality_factor, fs):
        b, a = iirnotch(notch_freq / (fs / 2), quality_factor)
        return filtfilt(b, a, data)
    
    data_init = np.copy(my_raw.EEG_data)
    time_to_pad = 1/notch_frequency*3
    n_to_pad = int(time_to_pad*my_raw.sampling_rate)
    
    for start, length in zip(my_raw.clean_runs, my_raw.clean_runs_lengths):
        data_to_filter = my_raw.EEG_data[start:start+length]
        data_padded = np.pad(data_to_filter, n_to_pad, mode='reflect')
        filtered_segment = iir_notch_filter(data=data_padded, notch_freq=notch_frequency, 
                                         quality_factor=q_factor, fs=my_raw.sampling_rate)
        filtered_segment = filtered_segment[n_to_pad:-n_to_pad]
        if len(data_to_filter) != len(filtered_segment):
            raise ValueError("The length of the EEG data segment and the filtered segment should be the same.")
        my_raw.EEG_data[start:start+length] = filtered_segment

    if np.array_equal(my_raw.EEG_data, data_init):
        raise ValueError("After notch filtering, there was no change to the data.")
    
    return my_raw

def preprocess(ar_obj):
    """
    Apply comprehensive preprocessing pipeline to generate multiple frequency bands.
    Implements clean-run detection and multi-band filtering to create frequency-specific
    versions of EEG data suitable for spectral analysis.
    
    Parameters:
    -----------
    ar_obj : alex_raw_efficient object
        Artifact-rejected EEG object with NaN-marked artifacts
        Must contain EEG_data array and sampling_rate
        
    Returns:
    --------
    list of alex_raw_efficient objects
        [notch_filtered, delta, theta, alpha, beta, trimmed_original]
        Each object contains filtered version of input data
        
    Processing Pipeline:
    -------------------
    
    Step 1: Clean Run Detection
    --------------------------
    Identifies continuous segments of valid (non-NaN) data:
    
    ```
    Signal:  [data, data, NaN, NaN, data, data, data, NaN]
    Runs:    [----run1----]         [-----run2-----]
    ```
    
    Algorithm:
    1. Create boolean mask: valid_data = ~np.isnan(EEG_data)
    2. Find transitions: diff(valid_data) identifies run boundaries
    3. Extract start/end indices of each continuous run
    4. Handle edge cases (runs at signal start/end)
    
    Step 2: Minimum Segment Filtering  
    ---------------------------------
    Converts segments shorter than min_seg_length to NaN:
    
    Rationale:
    - Short segments create filter artifacts at boundaries
    - Insufficient data for reliable frequency analysis
    - Butterworth filters need settling time
    
    Threshold: min_seg_length seconds (configurable)
    
    Step 3: Clean Run Recalculation
    -------------------------------
    After removing short segments, recalculates:
    - clean_runs: Array of starting indices for valid segments
    - clean_runs_lengths: Array of segment lengths (in samples)
    
    These attributes guide subsequent filtering operations.
    
    Step 4: Multi-Band Filtering
    ----------------------------
    Creates frequency-specific representations:
    
    1. Notch Filtered: Full spectrum with power line removal
       - Removes 50/60 Hz interference and extracts passband
    2. Delta Band (0.5-4 Hz):
    3. Theta Band (4-8 Hz):
    4. Alpha Band (8-13 Hz):
    5. Beta Band (13-30 Hz):
    6. Trimmed Original:
       - Unfiltered but segmented
       - Reference for filter validation
       - Preserves original spectral content
    
    Quality Assurance:
    -----------------
    Each filtered object maintains:
    - Original temporal structure (NaN positions preserved)
    - Clean run boundaries (for validation)
    - Processing history (change_log updates)
    - Metadata consistency (sampling rates, timestamps)
    
    Configuration Dependencies:
    --------------------------
    - min_seg_length: Minimum segment duration for analysis
    - full_range: Broadband frequency limits
    - delta_range, theta_range, alpha_range, beta_range: Band-specific limits
    - bandpass_order: Filter steepness
    - notch_frequency, notch_q: Power line removal parameters
    """
    window_size = vars_dict["min_seg_length"]

    # Step 1: Initialize important variables
    EEG_data = ar_obj.EEG_data

    # Step 2: Clean Run Detection Algorithm
    # ====================================
    # Identify continuous segments of valid (non-NaN) data for filtering
    
    # Create boolean mask: True for valid data, False for NaN artifacts
    valid_data = ~np.isnan(EEG_data)
    # Detect transitions between valid and invalid data
    # diff() identifies boundaries: 1 = NaN→valid, -1 = valid→NaN
    transitions = np.diff(valid_data.astype(int))

    # Extract run boundaries from transition points
    run_starts = np.where(transitions == 1)[0] + 1
    run_ends = np.where(transitions == -1)[0] + 1
    # Handle edge cases: runs that start/end at signal boundaries
    if valid_data[0]: # Signal starts with valid data
        run_starts = np.insert(run_starts, 0, 0)
    if valid_data[-1]: # Signal ends with valid data  
        run_ends = np.append(run_ends, len(valid_data))
    # Calculate run lengths for minimum segment filtering
    run_lengths = run_ends - run_starts
    # Convert short runs to NaN
    for start, length in zip(run_starts, run_lengths):
        if length < window_size*ar_obj.sampling_rate:
            EEG_data[start:start+length] = np.nan
    # Recalculate runs after converting short segments to NaN
    valid_data = ~np.isnan(EEG_data)
    transitions = np.diff(valid_data.astype(int))
    run_starts = np.where(transitions == 1)[0] + 1
    run_ends = np.where(transitions == -1)[0] + 1
    # Handle edge cases again
    if valid_data[0]:
        run_starts = np.insert(run_starts, 0, 0)
    if valid_data[-1]:
        run_ends = np.append(run_ends, len(valid_data))
        
    # Calculate final run lengths
    clean_runs = run_starts
    clean_runs_lengths = run_ends - run_starts
    # Add attributes to ar_obj
    ar_obj.EEG_data = EEG_data
    ar_obj.clean_runs = clean_runs
    ar_obj.clean_runs_lengths = clean_runs_lengths
    # Apply filters
    bandpassed_obj = bandpass_filter_signal(deepcopy(ar_obj), vars_dict["full_range"])
    notch_obj = notch_filter_signal(deepcopy(bandpassed_obj))
    delta_obj = bandpass_filter_signal(deepcopy(ar_obj), vars_dict["delta_range"])
    theta_obj = bandpass_filter_signal(deepcopy(ar_obj), vars_dict["theta_range"])
    alpha_obj = bandpass_filter_signal(deepcopy(ar_obj), vars_dict["alpha_range"])
    beta_obj = bandpass_filter_signal(deepcopy(ar_obj), vars_dict["beta_range"])
    return [notch_obj, delta_obj, theta_obj, alpha_obj, beta_obj, ar_obj]
     
def iterate_preprocess(save_base_dir):
    """
    Function that iterates through all _ar files matching specific number patterns and applies preprocessing steps
    """
    _ar_dir = os.path.join(os.path.split(save_base_dir)[0], r"2_artefact_rejected", r"artefact_rejected_files")

    # Clear directories if requested
    directories_to_clear = [
        save_base_dir,
        save_base_dir + r"\preprocessed_files_delta",
        save_base_dir + r"\preprocessed_files_theta",
        save_base_dir + r"\preprocessed_files_alpha",
        save_base_dir + r"\preprocessed_files_beta",
        save_base_dir + r"\preprocessed_files_original",
        save_base_dir + r"\preprocessed_files"
    ]
    
    for directory in directories_to_clear:
        os.makedirs(directory, exist_ok=True)
    
    for preprocessed_file_dir in directories_to_clear:
        if not os.path.exists(preprocessed_file_dir):
            continue
        for filename in os.listdir(preprocessed_file_dir):
            file_path = os.path.join(preprocessed_file_dir, filename)
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"Failed to delete {file_path}. Reason: {e}")



    # Iterate over artefact rejected files
    for _ar_file in os.listdir(_ar_dir):
        print(f"Processing {_ar_file}...")
        if _ar_file.endswith('.pkl'):
            import_path = os.path.join(_ar_dir, _ar_file)
            _ar_obj = import_pkl(full_file_path=import_path, parts_or_path="path")
            processed_arr = preprocess(_ar_obj)
            
            if processed_arr is not None:
                # Save processed objects (rest of the code remains the same)
                notch_obj = processed_arr[0]
                notch_obj.save_as_easy(save_base_dir + r"\preprocessed_files", '_pp')
                del notch_obj
                
                delta_obj = processed_arr[1]
                delta_obj.save_as_easy(save_base_dir + r'\preprocessed_files_delta', "_ppd")
                del delta_obj
                
                theta_obj = processed_arr[2]
                theta_obj.save_as_easy(save_base_dir + r'\preprocessed_files_theta', "_ppt")
                del theta_obj
                
                alpha_obj = processed_arr[3]
                alpha_obj.save_as_easy(save_base_dir + r'\preprocessed_files_alpha', "_ppa")
                del alpha_obj
                
                beta_obj = processed_arr[4]
                beta_obj.save_as_easy(save_base_dir + r'\preprocessed_files_beta', "_ppb")
                del beta_obj
                
                trimmed_obj = processed_arr[5]
                trimmed_obj.save_as_easy(save_base_dir + r'\preprocessed_files_original', "_ppo")
                del trimmed_obj
            else: #Case where all data was artefact
                continue