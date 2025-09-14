import os
import sys
current_dir = os.path.dirname(os.path.abspath(__file__))
utils_dir = os.path.join(current_dir, '..', 'data',"nicolet_files")
sys.path.append(utils_dir)


vars_dict = {    
    #File structure variables; recommended not to touch unless one is comfortable with programming
    "nicolet_folder":os.path.join(current_dir, '..', 'data',"nicolet_files"),
    "OBM_folder":os.path.join(current_dir, '..', 'data',"obm_files"),
    "base_dir": os.path.join(current_dir, '..', 'processing_run'),
    "nicolet_recording_time_data_path": os.path.join(current_dir, '..', 'data',"patient_data",'nicolet_file_times.xlsx'),
    "OBM_recording_time_data_path":os.path.join(current_dir, '..', 'data',"patient_data",'obm_file_times.xlsx'),
    "surgery_times_path": os.path.join(current_dir, '..', 'data',"patient_data",'surgery_time_data.xlsx'),
    "patient_data_path": os.path.join(current_dir, '..', 'data',"patient_data",'patient_data.xlsx'),
    "clear": True, #If true, we will delete the fis_obj dictionary used when building stats_gui

    #Basics of processing
    "windows":[[44,48]], #[[a,b]] where a is the start time of the postop window in hours postop, and b is the end time of the postop window in hours postop
    "nicolet_channels":['EEG P3-REF', 'EEG P4-REF'],
    "OBM_channel":["CrossEeg"],
    "downsample_rate":200,

    #Frequency bands
    "full_range":[0.5,30],
    "delta_range":[0.5,4],
    "theta_range":[4,7],
    "alpha_range":[7,13],
    "beta_range":[13,30],


    #Artefact rejection parameters
    "proportion_cutoff":1, #Proportion of signal that is allowed to be artefact; signals with more artefact than this are not processed further
    "time_cutoff":0, #1-proportion cutoff
    "zero_run_threshold":1, #Maximum number of SECONDS we are allowed to have a 0 run for before marking as artefact (in NEURAL_PY this was 1)
    "high_amp_collar":10, #Time in SECONDS for high-amplitude artefact collar (10 in NEURAL_PY)
    "high_amp_threshold":1000, #Voltages in uV above or below which is considered artefact (1500 in NEURAL_PY)
    "continuous_collar":0.5, #Time in SECONDS for continuous value artefact collar (0.5 in NEURAL_PY)
    "continuous_threshold":0.1, #Maximum time in SECONDS allowed for a continuous value before it is considered artefact (0.1 in NEURAL_PY)
    "sj_collar":0.5, #Time in SECONDS for sudden-jump artefact collar (0.5 in NEURAL_PY)
    "sj_threshold":200, #Maximum voltage difference in uV allowed for two consecutive values before they are considered artefact (200 in NEURAL_PY)

    #Preprocessing parameters
    "min_seg_length":20, #Minimum length in SECONDS that is required of a clean run to be kept; segments with a length under this are marked entirely as artefact before preprocessing
    "notch_frequency":50, #Frequency at which the notch filter should be applied. 50 Hz in Europe, usually 60 Hz in the US
    "notch_q":50, #q-factor for notch filter (higher values provide a tighter filter around notch_frequency)
    "bandpass_order":4 #Filter order for the passband filter used
}