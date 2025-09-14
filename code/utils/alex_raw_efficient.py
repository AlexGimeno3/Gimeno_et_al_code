"""
This is a variation of the alex_raw class that is built to contain everything the original one does, but without copying over everything from the MNE class (which I imagine is taking up a lot of space). It will also have other efficiency improvements.

Of note: this class is built to take an alex_raw object that has already been stitched (essentially, this is a class that is built to work after the best hour has been determined).
"""
import os
from datetime import datetime
import pickle
from WindowBuildError import WindowBuildError
import numpy as np
import math
import re
import gc

class alex_raw_efficient:
    def __init__(self, original_raw):
        """
        Initializes our class by taking an mne raw object (this way, we inherit all the mne_raw methods and attrbutes, and in theory we should still be able to use MNE documentation for the raw object). NB: EEGtimes_from_zero are reported as an array of milliseconds.
        
        Inputs:
        - mne_raw (tuple): A raw MNE object (the one we will be analyzing and working with)
        - viewTime1 (int): This is for visualization only. It is the time, IN SECONDS, that we want our viewing outputs to start
        - viewTime2 (int): This is for visualization only. It is the time, IN SECONDS, that we want our viewing outputs to end
        - which_channel (arr of str): This contains the channel of analysis we are interested in. If there are two channels here, we will create a bipolar signal of which_channel[0]-which_channel[1]
        - filePath (str): the file from which we loaded the EDF file
        - call_from (str): if the class is being instantiated from a certain function, this will be updated to reflect that specific function's needs. Only option is currently "postop_file_finder", which alters the channel choosing behavior to allow for easier automation.
        """
        
        #General info
        self.info = original_raw.info
        self.window_type = None
        if hasattr(original_raw, 'window_type'):
            self.window_type == original_raw.window_type
        else:
            self.window_type == "postop"
        

        #Signal information
        self.EEG_data = np.round(original_raw.EEG_data,3)
        self.EEG_times_from_zero = np.round(original_raw.EEG_times_from_zero,3)
        self.EEG_times_from_analyze_end_time = np.round(original_raw.EEG_times_from_analyze_end_time,3)
        self.sampling_rate = original_raw.sampling_rate
        self.nyquist_frequency = self.sampling_rate/2
        self.online_high_pass = original_raw.online_high_pass
        self.online_low_pass = original_raw.online_low_pass
        self.windows = original_raw.windows
        self.chunk = None
 
        #File storage data
        self.file_path_original = original_raw.file_path_original
        self.file_path_current = os.path.join(r"E:\Final_EEG_Pipeline\stitched_cut\stitched_cut_efficient", f"{original_raw.file_name_whole}.pkl")
        self.file_dir_current = r"E:\Final_EEG_Pipeline\stitched_cut\stitched_cut_efficient"
        self.saved_once = True

        #Other misc
        self.change_log = original_raw.change_log
        self.channels = original_raw.channels
        self.rec_start_date_time = original_raw.rec_start_date_time

    @property
    def signal_length_secs(self):
        """Recalculate and return the signal length in SECONDS."""
        return len(self.EEG_data)/self.sampling_rate

    @property
    def file_name_whole(self):
        my_str = os.path.splitext(os.path.basename(self.file_path_current))[0] #will return the entire file name, just without .pkl (e.g., FIS101post_P3-P4_stitched_cut_artefact_rejected_20_24_hrs)
        return my_str
    
    @property
    def file_name(self):
        return self.file_name_whole
    
    @property
    def analyze_end_time(self):
        return self.rec_start_date_time
        
    @property
    def num_fis(self):
        def fis(self):
            base_str = self.file_path_original
            match = re.search(r"(FIS\d+)", base_str)  # Use re.search to find the first match
            if match:
                fis_part = match.group(1)  # Extract the matched part
                return fis_part
            else:
                raise ValueError("This file's FIS could not be extracted from its file name.")
        
        base_str = fis(self)
        return ''.join(char for char in base_str if char.isdigit())  # Extract digits (e.g., 102 from FIS102)


    @property
    def fis_num(self):
        return int(self.num_fis)  # Extract digits (e.g., 102 from FIS102)
    
    @property
    def fis(self): #Redundant in case of confusion
        return self.num_fis
    
    @property
    def chunks(self):
        return self.windows
    
    @property
    def channel_name(self): #This is used for displaying names in graphs, etc
        if len(self.channels) == 1:
            parts = self.channels[0].split()  # Splits into ["EEG", "F3-REF"]
            electrode = parts[1].replace("-REF", "-reference")  # Converts "F3-REF" to "F3-reference"
            return electrode
        elif len(self.channels) == 2:
            str1 = self.channels[0]
            str2 = self.channels[1]
            electrode1 = str1.split()[1].replace("-REF", "")
            electrode2 = str2.split()[1].replace("-REF", "")
            # Combine the electrodes in "AB-CD" format
            return f"{electrode1}-{electrode2}"
    
    @property
    def num_windows(self):
        return len(self.windows)
   
    def add_change_log(self, changeText, suppress=False):
        """
        This method is used to write to an object's change-log every time something is done. This method must be called after each modification we want to log.

        Inputs:
        - changeText: the text (as a string) that will be added into the changelog

        Outputs: none (this method only updates our changelog)
        """
        now = datetime.now()
        datetime_str = now.strftime("%Y-%m-%d %H:%M:%S")
        if suppress == False:
            print(datetime_str+"\n"+changeText)
        self.change_log = self.change_log + "\n\n"+datetime_str+"\n"+changeText

    def save_as_easy(self, directory, suffix):
        """
        Function that helps us save a file while adding suffix. self.save_as has been optimized to work with the EDF reading step, whereas this function is good for everything else.

        Inputs:
        - directory (str) = the path where we would like the file to be dropped 
        - suffix (str): the suffix we would like to add onto the end of the file name (but before the .pkl file ending)
        """
        # Create directory if it doesn't exist
        if not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)
            self.add_change_log(f"Created directory {directory}", suppress=True)
        
        current_path = self.file_path_current
        file_name_no_pkl = os.path.splitext(os.path.basename(current_path))[0]
        new_file_name = file_name_no_pkl+suffix+'.pkl'
        new_path = os.path.join(directory, new_file_name)
        self.file_path_current = new_path
        self.add_change_log(f"File saved to {new_path}", suppress=True)
        with open(new_path, 'wb') as file:
            pickle.dump(self, file, protocol=pickle.HIGHEST_PROTOCOL)
        gc.collect()
        print(f"File {self.file_path_current} saved.")
    
    def time_to_index(self, time, times_arr=None, round_option=None, round_threshold=None):
        """
        Function that takes a given time value and returns the index associated with that data point in the chosen times_arr.

        Inputs:
        - time (float): time (i.e., recording start) IN SECONDS whose index we want to find.
        - times_arr (float arr): array of times IN SECONDS wherein we want to find time
        - round_option (str): can be None, "nearest", "up", or "down". If None, in cases where there is no exact match for time in time_from_zeros_array, throw an error.
        - round_threshold (float): the maximum amount of time in SECONDS allowed for time to differ from the time at the returned index.

        Outputs:
        - index (int): the index at which time can be found in self.EEG_times_from_zero.
        """


        round_error = -math.floor(math.log10(1/self.sampling_rate))  # Find number of decimal places
        
        # Default round_threshold to 1/sampling_rate if not provided
        if times_arr is None:
            times_arr = self.EEG_times_from_zero
        if round_threshold is None:
            round_threshold = 1 / self.sampling_rate
        eeg_times = np.array(times_arr)
        differences = np.abs(eeg_times - time)
        if round_option is None:
            if time in eeg_times:
                return np.where(eeg_times == time)[0][0]  # Return index of exact match
            else:
                raise ValueError(f"Error: No exact match found for time {time} in EEG_times_from_zero.")
        
        # Rounding logic
        elif round_option == "nearest":
            min_diff = round(np.min(differences),round_error)
            if min_diff <= round_threshold:
                return np.argmin(differences)  # Index of minimum difference
            else:
                print(f"Error: No time within {round_threshold * 1000} ms of {time}.")
                return None

        elif round_option == "up": #NB: assumes self.EEG_times_from_zero is ASCENDING (a fair assumption in our data)
            indices_up = np.where(eeg_times >= time)[0]
            if len(indices_up) > 0:
                index = indices_up[0]
                if (round(eeg_times[index] - time,round_error)) <= round_threshold:
                    return index
            print(f"Error: No time found within {round_threshold * 1000} ms after {time}.")
            return None
        
        elif round_option == "down":
            indices_down = np.where(eeg_times <= time)[0]
            if len(indices_down) > 0:
                index = indices_down[-1]
                if (round(time - eeg_times[index],round_error)) <= round_threshold:
                    return index
            print(f"Error: No time found within {round_threshold * 1000} ms before {time} ms.")
            return None
        else:
            print("Error: Invalid rounding option. Must be one of None, 'nearest', 'up', or 'down'.")
            return None
    
    def get_window_indices(self, postop_windows, times_from_postop = None):
        """
        This function is built to return indices in the self.EEG_times_from_analyze_end_time array corresponding to the postoperative windows of interest.

        Inputs:
        - postop_windows (arr of int or float arr): this is an array of the time windows we want in HOURS. For example, [[20, 24],[44, 48]] represents those time windows at 20hrs-24hrs postop as well as those times at 44-48 hours postop.
        - times_from_postop_arr: set to self.EEG_times_from_analyze_end_time by default; this is built in the functions in stitch_nicolet_files.py

        Outputs:
        - indices_arr (arr of int arr): this is an array of indices corresponding to our ranges. So [[6000, 6500],[7000, 7600]] indicates that the first time window is contained in the data at self.EEG_times_from_analyze_end_time[6000] to self.EEG_times_from_analyze_end_time[6500], and that the second window's data is contained at similar time points.
        - windows_arr (arr of int arr): this is the array of windows that the EEG_data at least partially overlapped with
        """
        windows_arr = []
        if times_from_postop == None:
            try:
                times_from_postop = self.EEG_times_from_analyze_end_time
            except AttributeError:
                print("There is no array documenting times postoperatively. Make sure this is a stitched Nicolet file (or a processed non-Nicolet file) and try again.")
                raise
        windows_in_seconds = [[n1 * 3600, n2 * 3600] for n1, n2 in postop_windows]

        indices_arr = [] #The indices here are indices in the self.EEG_times_from_analyze_end_time. An index of 0 refers to the beginning of the recording; and index of len(times_from_postop)-1 refers to the end of the recording.
        i = 0
        for window in windows_in_seconds:
            time1 = window[0] #time1 and time2 refer to times from postop, NOT times from zero
            time2 = window[1]
            rec_time_1 = times_from_postop[0]
            rec_time_2 = times_from_postop[-1]
            if (time1<=rec_time_1 and time2 <=rec_time_1) or (time1>=rec_time_2 and time2 >=rec_time_2): #This takes care of our cases where our window is not in the recording at all
                i = i+1
                continue #Will not add anything to add_arr
            elif rec_time_1<=time1 and rec_time_2>=time2: #i.e., window falls entirely within our recording
                index1 = self.time_to_index(time1,times_arr=times_from_postop,round_option="up")
                index2 = self.time_to_index(time2,times_arr=times_from_postop,round_option="down")
            elif rec_time_1>=time1 and rec_time_2<=time2: #i.e., our window is larger than our recording; i.e., the entire recording is in the window
                #Return the indices of the whole
                index1 = 0 
                index2 = len(times_from_postop)-1
            elif rec_time_1<=time1 and rec_time_2<=time2: #i.e., the beginning of our window falls within the recording, but the end of it falls outside
                index1 = self.time_to_index(time1,times_arr=times_from_postop,round_option="up")
                index2 = len(times_from_postop)-1 #to get last index of times_from_postop
            elif rec_time_1>=time1 and rec_time_2>=time2: #i.e., the beginning of our window falls outside the recording, but the end of it falls within the recording
                index1 = 0
                index2 = self.time_to_index(time2,times_arr=times_from_postop,round_option="down") #to get last index of times_from_postop
            else:
                i = i+1
                continue #Will not add anything to add_arr
            add_arr = [index1, index2]
            indices_arr.append(add_arr)
            windows_arr.append(postop_windows[i])
            i = i+1
        if len(indices_arr)==0:
            raise WindowBuildError(
                f"None of this file was in the windows of interest (hr): "
                f"{[[(elem / 3600.0) for elem in win] if isinstance(win, (list, tuple)) else win / 3600.0 for win in windows_in_seconds]}. "
                f"The file started {round(self.EEG_times_from_analyze_end_time[0]/3600,2)} hours after the end of surgery "
                f"(which was on {self.analyze_end_date_time.strftime('%d %B, %Y %H:%M')}) "
                f"and ran {round(self.signal_length_secs/3600,2)} hrs. "
                f"This often happens, for example, when we look at a surgery file, but we are interested in postop files."
            )
        else:
            return indices_arr, windows_arr