import pandas as pd
from datetime import time
from stat_vars import test_dict, wmi_dict, bis_dict, sex_dict, dx_dict
import numpy as np
import math
from tqdm import tqdm
import gc
import ast
import os
import sys
current_dir = os.path.dirname(os.path.abspath(__file__))
utils_dir = os.path.join(current_dir, '..', '..', 'utils')
sys.path.append(utils_dir)
from utils import config
vars_dict = config.get_vars_dict()

class fis_object():
    def __init__(self, fis, ts_dict, EEG_data, EEG_times, op_status, srate):
        """
        Initializes our class, whose goal is to efficiently store and return time series data at different time points. There will be one of these per FIS.

        Inputs:
        - fis (int): the fis of this stats object
        - ts_dict (dict): a dictionary of var_objects for this given fis

        Useful attributes:
        - filters: self.ND_filter, self.ecc_filter, self.machine_type_filter, self.ignore_filter
        - self.ECC_times: {"ecc_start": None, "clamp_start": None, "clamp_end": None, "ecc_end": None} 
        """
        
        #Step 1: initialize filters
        self.fis = fis
        self.ts_dict = ts_dict        
        self.zero_point = "end"
        self.ND_filter = 1
        self.ecc_filter = 1
        self.machine_type_filter = 1
        self.ignore_filter = 1
        self.artefact_filter = 1
        self.seizure_filter = 1
        self.cx_start_time = None
        self.cx_end_time = None
        for var_obj in list(ts_dict.values()):
            # Check for repeats in the times array
            times = var_obj.times
            if len(times) != len(set(times)):
                raise ValueError(f"Times should never repeat, but this happened in FIS{self.fis}, variable {var_obj.name}.")
        self.EEG_data = EEG_data
        self.EEG_times = EEG_times
        self.skip_str = None #When not None, indicates the reason the file was skipped in analysis
        self.analyzed_whole_sample = 0 #When 1, indicates that this FIS was indeed analyzed in the whole sample
        self.analyzed_seizure_sample = 0
        self.analyzed_across_seizures = 0
        self.op_status = op_status
        self.sampling_rate = srate
        
        
        #Step 2: Get ECC data
        #Step 2.1: Get ECC and clamp start/end times
        self.cpb_present = 0
        self.ECC_times = {"ecc_start": None, "clamp_start": None, "clamp_end": None, "ecc_end": None}
        ecc_excel_dir = vars_dict["patient_data_path"]
        df = pd.read_excel(ecc_excel_dir)
        # Check if the patient with fis_num equal to fis has cpb_used equal to 1
        if df.loc[df['fis_num'] == fis, 'cpb_used'].values[0] == 1:
            self.cpb_present = 1
            # Open the surgery times Excel sheet
            excel_path = vars_dict["surgery_times_path"]
            df = pd.read_excel(excel_path)
            # Find the row with matching fis_num
            row = df[df['fis_num'] == int(fis)]
            check_strs = ["time_start_ecc", "time_start_clamp","time_end_clamp", "time_end_ecc"]
            for col_name, dict_key in zip(check_strs, self.ECC_times.keys()):
                time_start_ecc_value = str(df.loc[df['fis_num'] == fis, col_name].values[0])
                # Check if the value is in the format XXXX (4 digits)
                if isinstance(time_start_ecc_value, str) and len(time_start_ecc_value) == 4:
                    # Convert to hours and minutes
                    hours = int(time_start_ecc_value[:2])
                    minutes = int(time_start_ecc_value[2:])
                    # Create a time object
                    ecc_start_time = time(hour=hours, minute=minutes)
                    # Store it in self.ECC_times dictionary
                    self.ECC_times[dict_key] = ecc_start_time
        else:
            pass
        #Step 2.2: Get surgery start/end times
        excel_path = vars_dict["surgery_times_path"]
        df_surgery = pd.read_excel(excel_path)
        surgery_start_value = f"{int(df_surgery.loc[df_surgery['fis_num'] == fis, 'time_start_cx'].values[0]):04d}"
        surgery_end_value = f"{int(df_surgery.loc[df_surgery['fis_num'] == fis, 'time_end_cx'].values[0]):04d}"
        if isinstance(surgery_start_value, str) and len(surgery_start_value) == 4 and surgery_start_value.isdigit():
            start_hours = int(surgery_start_value[:2])
            start_minutes = int(surgery_start_value[2:])
            self.cx_start_time = time(hour=start_hours, minute=start_minutes)
        else:
            raise ValueError("No surgery start time found.")
        if isinstance(surgery_end_value, str) and len(surgery_end_value) == 4 and surgery_end_value.isdigit():
            end_hours = int(surgery_end_value[:2])
            end_minutes = int(surgery_end_value[2:])
            self.cx_end_time = time(hour=end_hours, minute=end_minutes)
        else:
            raise ValueError("No surgery end time found.")
        #Step 2.3: Convert ECC times to times relative to surgery end in SECONDS
        for key in self.ECC_times:
            if isinstance(self.ECC_times[key], time):
                # Calculate the time difference in hours
                ecc_time = self.ECC_times[key]
                end_time = self.cx_end_time
                # Convert both times to minutes since midnight
                ecc_minutes = ecc_time.hour * 60 + ecc_time.minute
                end_minutes = end_time.hour * 60 + end_time.minute
                # Handle cases where ECC time might be on the next day
                if ecc_minutes < end_minutes:
                    # Assume both times are on the same day
                    diff_minutes = ecc_minutes - end_minutes
                else:
                    # If ECC time is later in the day than end time, assume it's from the same day
                    diff_minutes = ecc_minutes - end_minutes
                # Convert difference to seconds
                diff_hours = diff_minutes*60
                # Store the difference in hours
                self.ECC_times[key] = diff_hours
        if isinstance(self.cx_start_time, time) and isinstance(self.cx_end_time, time):
            start_minutes = self.cx_start_time.hour * 60 + self.cx_start_time.minute
            end_minutes = self.cx_end_time.hour * 60 + self.cx_end_time.minute
            diff_minutes = start_minutes - end_minutes  # Will be negative if start is before end
            self.cx_start_time = diff_minutes*60  # Convert to secs
        self.cx_end_time = 0

        #Step 3: Get ND status
        # Open the Excel sheet to get the Final_ND value
        eeg_data_path = vars_dict["patient_data_path"]
        df_eeg = pd.read_excel(eeg_data_path)
        # Find the row in "fis_num" that matches self.fis
        matching_row = df_eeg[df_eeg['fis_num'] == self.fis]
        if not matching_row.empty:
            # Get the value in the "Final_ND" column
            final_nd_value = matching_row['Final_ND'].values[0]
            # Set self.ND based on the value
            if final_nd_value == 0:
                self.ND = 0
            elif final_nd_value == 1:
                self.ND = 1
            else:
                self.ND = None
        else:
            # If no matching row is found
            self.ND = None
    
        
        def import_clean(matching_row, col_name, data_type):
            """
            Function that imports and returns the item in matching_row and col_name of the given data_type. Returns the converted value, or None if value is "#NULL!" or NaN. Raises ValueError if conversion to the specified data_type fails.
            """
            # Get the value from the matching row
            value = matching_row[col_name].item()
            # Check for null values
            if value == "#NULL!" or pd.isna(value):
                return None
            # Attempt conversion
            try:
                return data_type(value)
            except (ValueError, TypeError) as e:
                raise ValueError(f"Failed to convert '{value}' to {data_type.__name__} for column '{col_name}'") from e
        
        #Step 4 : Get machine type
        eeg_data_path = vars_dict["patient_data_path"]
        df_eeg = pd.read_excel(eeg_data_path)
        # Find the row in "fis_num" that matches self.fis
        matching_row = df_eeg[df_eeg['fis_num'] == self.fis]
        if not matching_row.empty:
            # Get the value in the "machine_type" column
            mt = int(matching_row['machine_type'].values[0])
            # Set self.ND based on the value
            if mt == 1:
                self.machine_type = "nicolet"
            elif mt == 2:
                self.machine_type = "obm"
            else:
                self.machine_type = None

            #Step 5: Import other important data
            #Step 5.1: Import demo/clinical data
            date_offset = (self.cx_start_time + self.EEG_times[0])/3600 #cx_start_time is in seconds from 0000
            date_offset = math.floor(date_offset/24)
            
            self.ga_eeg = import_clean(matching_row, "gestational_age_at_surgery", int)
            if not self.ga_eeg is None:
                self.ga_eeg = self.ga_eeg + date_offset
            self.pma_eeg = import_clean(matching_row, "dol_at_surgery", int)
            if not self.pma_eeg is None:
                self.pma_eeg = self.pma_eeg + date_offset
            self.sex = import_clean(matching_row, "sex", int)
            self.maternal_age = import_clean(matching_row, "maternal_age_at_birth", int)
            self.birth_weight = import_clean(matching_row, "weight_at_birth", int)
            self.cardiac_dx = import_clean(matching_row, "cardiac_diagnosis", int)
            self.sts_eacts_category = import_clean(matching_row, "STS_EACTS_category", int)
            self.aristotle_score = import_clean(matching_row, "Aristoteles_score", int)
            self.aristotle_category = import_clean(matching_row, "Aristoteles_category", int)
            self.rachs1_category = import_clean(matching_row, "RACHS1_category", int)
            #Step 5.1.1: Import CPB data
            self.cpb_used = import_clean(matching_row, "cpb_used", int)
            self.cpb_time = import_clean(matching_row, "cpb_time", int)
            self.cpb_clamp_used = import_clean(matching_row, "cpb_clamp_used", int)
            self.cpb_clamp_time = import_clean(matching_row, "cpb_clamp_time", int)
            self.circulatory_arrest_used =import_clean(matching_row, "circulatory_arrest_used", int)
            self.circulatory_arrest_time = import_clean(matching_row, "circulatory_arrest_time", int)
            self.clamp_used_coarct_repair = import_clean(matching_row, "aortic_clamp_used_coarct_repair", int)
            self.clamp_time_coarct_repair = import_clean(matching_row, "aortic_clamp_time_coarct_repair", int)
            self.uf_used = import_clean(matching_row, "ultrafiltration_used", int)
            self.uf_volume = import_clean(matching_row, "ultrafiltration_volume", int)
            self.acp_used = import_clean(matching_row, "acp_used", int)
            self.acp_flow_rate = import_clean(matching_row, "acp_flow_rate", int)
            self.acp_time = import_clean(matching_row, "acp_time", int)
            self.cooling_temperature = import_clean(matching_row, "cooling_temperature", int)
            self.hypothermia_time = import_clean(matching_row, "hypothermia_time", int)
            #Step 5.2: Import test-related data
            self.test_used = import_clean(matching_row, "test_used", int)
            self.dol_for_test_used = import_clean(matching_row, "DoL_for_test_used", int)
            #Step 5.2.1: Bayley stats
            self.bayley_cognitive_composite = import_clean(matching_row, "bayley_cognitive_composite", int)
            self.bayley_language_composite = import_clean(matching_row, "bayley_language_composite", int)
            self.bayley_motor_composite = import_clean(matching_row, "bayley_motor_composite", int)
            #Step 5.2.2: VL stats
            self.vl_abc = import_clean(matching_row, "vineland_adaptive_behavior_composite", int)
            self.vl_communication = import_clean(matching_row, "vineland_communication_total", int)
            self.vl_daily_activities = import_clean(matching_row, "vineland_daily_activities_total", int)
            self.vl_socialization = import_clean(matching_row, "vineland_socialization_total", int)
            self.vl_psychomotor_development = import_clean(matching_row, "vineland_psychomotor_development_total", int)
            #Step 5.3: Import NIRS data
            self.time_nirs_under_40 = import_clean(matching_row, "time_under_40_percent", int)
            self.prop_nirs_under_40 = import_clean(matching_row, "prop_under_40_percent", int)
            self.time_nirs_under_50 = import_clean(matching_row, "time_under_50_percent", int)
            self.prop_nirs_under_50 = import_clean(matching_row, "prop_under_50_percent", int)
            self.time_nirs_out_of_50_70 = import_clean(matching_row, "time_out_of_50_70_percent", int)
            self.prop_nirs_out_of_50_70 = import_clean(matching_row, "prop_out_of_50_70_percent", int)
            self.time_nirs_over_85 = import_clean(matching_row, "time_over_85_percent", int)
            self.prop_nirs_over_85 = import_clean(matching_row, "prop_over_85_percent", int)
            #Step 5.4: Import biomarker data
            self.preop_enolase = import_clean(matching_row, "preop_enolase", float)
            self.preop_S100B = import_clean(matching_row, "preop_S100B", float)
            self.postop_enolase = import_clean(matching_row, "postop_enolase", float)
            self.postop_S100B = import_clean(matching_row, "postop_S100B", float)
            self.enolase_24_hr = import_clean(matching_row, "@24_hr_enolase", float)
            self.S100B_24_hr = import_clean(matching_row, "@24_hr_S100B", float)
            self.enolase_72_hr = import_clean(matching_row, "@72_hr_enolase", float)
            self.S100B_72_hr = import_clean(matching_row, "@72_hr_S100B", float)
            #Step 5.5: Import MRI data
            self.stroke = import_clean(matching_row, "stroke", int)
            self.wmi = import_clean(matching_row, "white_matter_injury", int)
            self.wmi_details = import_clean(matching_row, "white_matter_injury_details", int)
            self.cerebellar_hemorrhage = import_clean(matching_row, "cerebellar_hemorrhage", int)
            self.subdural_hemorrhage = import_clean(matching_row, "subdural_hemorrhage", int)
            self.intraventricular_hemorrhage = import_clean(matching_row, "intraventricular_hemorrhage", int)
            self.global_hypoxic_lesion = import_clean(matching_row, "global_hypoxic_lesion", int)
            self.bis = import_clean(matching_row, "brain_injury_severity_score", int)
            #Step 5.6: import seizure data
            self.seizures_intraop = import_clean(matching_row, "seizures_intraop", int)
            self.seizures_postop = import_clean(matching_row, "seizures_postop", int)
        else:
            raise ValueError(f"File {self.fis} has an EEG recording but no reported machine type.")
        
        #--------------------------------------------- New stuff starts below here
        chunks_df = fis_object._load_chunk_df(self)
        # if self.op_status == "intraop":
        #     raise ValueError("Re-run the all code!")
        chunk_val = chunks_df.loc[((chunks_df["fis"] == str(self.fis)) | (chunks_df["fis"] == self.fis)) & (chunks_df["chunk"] != "-"), "chunk"].iloc[0]
        chunk_hrs = np.asarray(ast.literal_eval(chunk_val), dtype=float)
        #chunk_hrs = [20,24]
        chunk_sec = np.asarray(chunk_hrs) * 3600.0
        start_sec, end_sec = chunk_sec
        # --- sampling-rate check ------------------------------------------------------
        if len(self.EEG_times) > 1:
            dt = np.diff(self.EEG_times)
            fs_max = 1.0/np.min(dt)
            fs_min = 1.0/np.max(dt)
            if not np.isclose(fs_max, self.sampling_rate, rtol=1e-3, atol=1e-3) and np.isclose(fs_min, self.sampling_rate, rtol=1e-3, atol=1e-3):
                raise ValueError(f"Detected max sampling rate {fs_max:.4f} Hz and min sampling rate {fs_min:.4f} Hz ≠ expected {self.sampling_rate:.4f} Hz")
            else:
                #print(f"Max sampling rate: {fs_max:.4f} Hz; min sampling rate: {fs_min:.4f}.")
                pass
        # --- padding sizes (in samples) ----------------------------------------------
        dt_expected = 1.0 / self.sampling_rate
        pad_before = int(round(max(self.EEG_times[0] - start_sec, 0) / dt_expected))
        pad_after = int(round(max(end_sec - self.EEG_times[-1], 0) / dt_expected))
        # --- create padding -----------------------------------------------------------
        data_pad_before = np.full(pad_before, np.nan)
        data_pad_after = np.full(pad_after, np.nan)
        self.EEG_data = np.concatenate((data_pad_before, self.EEG_data, data_pad_after))
        # --- rebuild corresponding time-vector ----------------------------------------
        total_samples = len(self.EEG_data)
        self.EEG_times = start_sec + np.arange(total_samples) * dt_expected
        print(f"Padded {pad_before} samples before and {pad_after} samples after ({pad_before + pad_after} total).")
        self.build_artefact_dict()

        

    from functools import lru_cache
    @staticmethod
    @lru_cache(maxsize=None)
    def _load_eeg_df():
        return pd.read_excel(vars_dict["patient_data_path"])
    @staticmethod
    @lru_cache(maxsize=None)
    def _load_surg_df():
        return pd.read_excel(vars_dict["surgery_times_path"])
    @staticmethod
    @lru_cache(maxsize=None)
    def _load_chunk_df(obj):
        mach_str = "OBM" if obj.machine_type=="obm" else "nicolet"
        if vars_dict["windows"] == [["surgery"]]:
            path = os.path.join(vars_dict["base_dir"],"intraop_files","1_stitched_cut",f"{mach_str}_raw_file_timestamps.xlsx")
        else:
            path = os.path.join(vars_dict["base_dir"],"postop_files","1_stitched_cut",f"{mach_str}_raw_file_timestamps.xlsx")
        return pd.read_excel(path)
    
    @property
    def filter(self):
        return self.ND_filter == 1 and self.ecc_filter == 1 and self.machine_type_filter == 1 and self.ignore_filter == 1 and self.seizure_filter == 1 and self.artefact_filter == 1
    
    @property
    def ecc_status(self):
        if self.cpb_present == 1 and any(value is not None for value in self.ECC_times.values()):
            return 1
        else:
            return 0
    
    @property
    def test_str(self):
        if not self.test_used is None:
            return test_dict[self.test_used]
        else:
            return None
    
    @property
    def wmi_details_str(self):
        if not self.wmi_details is None:
            return wmi_dict[self.wmi_details]
        else:
            return None

    @property
    def bis_str(self):
        if not self.bis is None:
            return bis_dict[self.bis]
        else:
            return None
    
    @property
    def sex_str(self):
        if not self.sex is None:
            return sex_dict[self.sex]
        else:
            return None

    @property
    def dx_str(self):
        if not self.dx is None:
            return dx_dict[self.dx]
        else:
            return None
    
    def build_artefact_dict(self):
        """
        Output:
        - art_dict (dict): builds a dictionary, where the keys are tuples of time values IN SECONDS and the values are 5-item float arrays with the following information: [proportion_art, num_arts, time_of_art, first_art, last_art]
        """
        art_dict = {}
        if self.op_status == "intraop":
            #Case for entire surgery; last = True in these cases because there is no iterative processing of artefacts
            start_time = self.cx_start_time
            end_time = self.cx_end_time
            art_dict[(start_time,end_time)] = self.build_artefact_ret(start_time,end_time, last = True)
            #Case for ecc
            start_time = self.ECC_times["ecc_start"]
            end_time = self.ECC_times["ecc_end"]
            if not (start_time is None or end_time is None):
                art_dict[(start_time,end_time)] = self.build_artefact_ret(start_time,end_time, last=True)
            else:
                art_dict[(start_time,end_time)] = None
            #Case for clamp
            start_time = self.ECC_times["clamp_start"]
            end_time = self.ECC_times["clamp_end"]
            if not (start_time is None or end_time is None):
                art_dict[(start_time,end_time)] = self.build_artefact_ret(start_time,end_time, last=True)
            else:
                art_dict[(start_time,end_time)] = None
            #Case for half
            start_time = self.ECC_times["ecc_start"]
            end_time = self.ECC_times["ecc_end"]
            if not (start_time is None or end_time is None):
                end_time_half = start_time + (end_time - start_time) / 2
                art_dict[(start_time,end_time_half)] = self.build_artefact_ret(start_time,end_time, last=True)
            else:
                end_time_half = None
                art_dict[(start_time,end_time_half)] = None
        elif self.op_status == "postop":
            first_hr = math.floor(self.EEG_times[0]/3600)
            latest_hr = math.ceil(self.EEG_times[-1]/3600) #I.e., if we go to 72.3 hrs (in secs), i iterates from 0 to 72. This means latest_hr will be 73 (since range(0,latest_hr) will go 0 to 72)
            for i in tqdm(range(first_hr,latest_hr), desc=f"Building artefact dictionary..."):
                #Handle special case where we are in the last hour
                if i == (latest_hr - 1):
                    start_time = i*3600
                    end_time = self.EEG_times[-1]
                    art_dict[(start_time,end_time)] = self.build_artefact_ret(start_time,end_time, last=True)
                else: #Handle normal cases
                    start_time = i*3600
                    end_time = (i+1)*3600
                    art_dict[(start_time, end_time)] = self.build_artefact_ret(start_time,end_time, last=False)
        self.art_dict = art_dict
        del self.EEG_data, self.EEG_times
        gc.collect()

    def build_artefact_ret(self, time1, time2, last=False):
        """
        Returns the artefact info to be used by build_artefact_dict
        """
        if not last:
            mask = (self.EEG_times >= time1) & (self.EEG_times < time2)
        else:
            mask = (self.EEG_times >= time1) & (self.EEG_times <= time2)
        # Return the filtered values
        art_data = self.EEG_data[mask]
        if art_data.size == 0: #case where the entire hour isn't present. NB: no checking for artefacts has occurred yet, so this returns 0 only when there is no data in the time range of interest.
            raise ValueError(f"The times requested ({time1}, {time2} secs) are outside the range of the preprocessed files: ({self.EEG_times[0]}, {self.EEG_times[-1]} secs)")
        # Proportion of artefact (NaNs)
        proportion_artefact = np.isnan(art_data).sum() / len(art_data)
        # Detect transitions from non-NaN to NaN and vice versa
        isnan = np.isnan(art_data)
        # Pad with False to detect edges properly
        padded = np.pad(isnan, (1, 1), constant_values=False)
        # XOR detects changes in the boolean array
        changes = padded[1:] != padded[:-1]
        # Each pair of changes (True to False or False to True) marks one artefact segment
        num_artefacts = np.count_nonzero(changes) // 2
        time_artefact = proportion_artefact * (time2 - time1) / 3600
        first_art = np.isnan(art_data[0])
        last_art = np.isnan(art_data[-1])
        return [proportion_artefact, num_artefacts, time_artefact, first_art, last_art]
    
    def build_artefact_info(self, time1, time2):
        """
        - time1 and time2 are in HOURS from surgery end. These values come from user input, and we are doing checking here to ensure that they are only exact hours.
        """
        if self.op_status == "postop":
            time1 = time1*3600
            time2 = time2*3600
        if self.op_status == "postop": #intraop case does not have this same issue of needing to ensure monotonically increasing keys, so we only check in the postop case        
            keys_arr = list(self.art_dict.keys())
            keys_arr = [item for sub_list in keys_arr for item in sub_list]
            if not all(keys_arr[i+1]>=keys_arr[i] for i in range(len(keys_arr)-1)):
                raise ValueError(f"keys_arr should be monotonically increasing, but it is not: list(self.art_dict.keys()) = {list(self.art_dict.keys())}.")
        
        art_key = (time1, time2)
        if self.op_status == "intraop":
            if not art_key in list(self.art_dict.keys()):
                raise ValueError(f"We wanted to see artefact data from {(time1, time2)}, but the only available times were {list(self.art_dict.keys())}.")
            art_arr = self.art_dict[art_key]
            prop_art = art_arr[0]
            num_art = art_arr[1]
            time_art = art_arr[2]
        elif self.op_status == "postop":
            #Error check
            start_times = [my_key[0] for my_key in list(self.art_dict.keys())]
            end_times = [my_key[1] for my_key in list(self.art_dict.keys())]
            if not(time1 in start_times and time2 in end_times):
                raise ValueError(f"You asked to check artefacts in the time range {(start_times, end_times)}. Unfortunately, the start times we have are {start_times} and the end times we have are {end_times}.")
            prop_arr = []#This will be an array of 1x2 arrays; each 1x2 array will contain the proportion of artefact as well as a weight for that proportion (which will always be 1 save for the last case)
            num_arr = [] #This will be an array of 1x3 arrays. The first entry is the number of NaN runs in this slice; the second is if the previous segment ended with an artefact, and the third is if the current segment began with an artefact
            time_art = 0 #This will be a simple iterative addition
            num_art = 0
            if not (time1//3600 == time1/3600 and time2//3600 == time2/3600):
                raise ValueError(f"It appears the provided time values were not whole hours: time1 = {time1} secs, {time1/3600}; time2 = {time2} sec, {time2/3600} hrs")

            t1_hrs = int(time1//3600)
            t2_hrs = int(time2//3600)
            #Calculate proportion of artefacts and time_art
            for i in range(0,int(t2_hrs-t1_hrs)):
                #Special last case
                if i+t1_hrs+1 == t2_hrs:
                    #t2_hrs would be the 2nd item in the key, which may not be an integer. Therefore, we will look for the key whose first value equals i+t1_hrs
                    keys_list = list(self.art_dict.keys())
                    fin_ind = [index for index,value in enumerate(keys_list) if value[0] == (i+t1_hrs)*3600]
                    if fin_ind:
                        my_key = keys_list[fin_ind[0]] #Since this is in the postop case only (where the first item in each key is by definition increasing hour-on-hour), the correct slice is always chosen
                    else: 
                        raise ValueError("No fin_ind found.")
                    t_diff_prop = (my_key[1]-my_key[0])/3600
                    prop_arr.append([self.art_dict[my_key][0], t_diff_prop])
                else: #All other values
                    my_key = ((t1_hrs+i)*3600,(t1_hrs+i+1)*3600)
                    if not math.isclose((my_key[1]-my_key[0])/3600, 1, rel_tol=0.001):
                        raise ValueError(f"(my_key[1]-my_key[0])/3600 should be 1, but ended up {(my_key[1]-my_key[0])/3600}")
                    prop_arr.append([self.art_dict[my_key][0], (my_key[1]-my_key[0])/3600]) #1 because non-last entries in the dictionary are always for the entire hour (and so weighted fully)
                time_art += self.art_dict[my_key][2]
            
            #Calculate number of artefacts
            #First, build num_arr
            prev_key = ()
            for i in range(t2_hrs-t1_hrs):
                temp_arr = [None, None, None]
                #Assign key
                if i+t1_hrs+1 == t2_hrs: #Special last case
                    keys_list = list(self.art_dict.keys())
                    fin_ind = [index for index,value in enumerate(keys_list) if value[0] == (i+t1_hrs)*3600]
                    if fin_ind:
                        my_key = keys_list[fin_ind[0]] #Since this is in the postop case only (where the first item in each key is by definition increasing hour-on-hour), the correct slice is always chosen
                    else: 
                        raise ValueError("No fin_ind found.")
                else: #All other values
                    my_key = ((t1_hrs+i)*3600,(t1_hrs+i+1)*3600)
                #Logic starts here
                temp_arr[0] = self.art_dict[my_key][1] #Number of NaN runs in this hour
                temp_arr[2] = self.art_dict[my_key][4] #True if last val is NaN
                if not (i == 0):
                    temp_arr[1] = self.art_dict[prev_key][4] #This updates temp_arr to contain if the previous hour-long segment ended in an artefact or not. For the first segment, we will simply default to None. NB: since first iteration is excluded, prev_key will always be defined.
                num_arr.append(temp_arr)
                prev_key = my_key
            #Quality check
            if len(num_arr) != t2_hrs-t1_hrs:
                raise ValueError(f"num_arr had len {len(num_arr)}, but there were {t2_hrs-t1_hrs} segments!")
            #Logic to calculate final num_art
            #Special first case: add number of NaN runs to num_art
            #Middle cases and last case: add number of NaN runs to num_arr; if the previous entry in num_arr ended with NaN AND the current entry starts with NaN, add (number of NaN runs - 1) to num_arr
            for i in range(len(num_arr)):
                if i == 0:
                    num_art += num_arr[i][0]
                else:
                    num_art += num_arr[i][0] - num_arr[i][1]*num_arr[i][2] #Will subtract 1 from number of runs if first val is NaN in this run and last val of previous run was NaN (this last run would have already been counted)
            prop_art = sum(entry[0]*entry[1] for entry in prop_arr)/sum(entry[1] for entry in prop_arr)

        self.proportion_artefact = prop_art
        self.num_artefacts = num_art
        self.time_artefact = time_art
        
    
    def ts(self, var_name, time1, time2):
        """
        Method to return the values and time stamps for a given variable of interest in a certain time window.
        """
        ts_vals = self.ts_dict[var_name].vals
        ts_times = self.ts_dict[var_name].times
        mask = (ts_times >= time1) & (ts_times <= time2)
        # Return the filtered values
        return ts_vals[mask], ts_times[mask]