"""
This class handles the loading of .edf files and extraction of potentially useful data. It is then later on converted into an alex_raw_efficient object to reduce memory usage and allow for quicker file loading/processing/saving.
"""

import matplotlib
matplotlib.use('Qt5Agg')
matplotlib.rcParams['text.usetex'] = False
import matplotlib.pyplot as plt
from NEURAL_py_EEG import preprocessing_EEG
import os
from datetime import datetime
import pickle
import sys
import tkinter as tk
from tkinter import filedialog
file_import_path = r"E:\Final_EEG_Pipeline\utils"
sys.path.append(file_import_path)
from WindowBuildError import WindowBuildError
import numpy as np
import pygetwindow as gw
import pyautogui
import math
from scipy.signal import welch
import re

class alex_raw:
    def __init__(self, mne_raw, viewTime1=0, viewTime2=5, which_channel=None, filePath = None, call_from = None):
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
        self.__dict__.update(mne_raw.__dict__)  # Copy all attributes from the mne_raw object
        self.info = mne_raw.info

        def choose_channel(self):
            # Display all channels with indices
            print("Available channels:")
            for i, channel in enumerate(self.channels):
                print(f"{i + 1}: {channel}")
            # Prompt user to choose a channel by number
            while True:
                try:
                    choice = int(input("Enter the number of the channel you want to choose: ")) - 1
                    if 0 <= choice < len(self.channels):
                        print(f"You chose: {self.channels[choice]}")
                        return self.channels[choice]
                    else:
                        print(f"Please enter a number between 1 and {len(self.channels)}.")
                except ValueError:
                    print("Invalid input. Please enter a valid number.")
        
        self.call_from = call_from
        self.file_path_original = filePath
        self.file_path_current = filePath #This includes everything, including the file extension
        self.file_dir_current = ""
        self.saved_once = False #Once the file has been saved once, this is turned to True. This allows us to manage the file name/extension in the save_as function

        self.change_log = ""
        #Add initialization to change_log
        now = datetime.now()
        datetime_str = now.strftime("%Y-%m-%d %H:%M:%S")
        self.change_log = datetime_str+"\n"+ f"File built on "+datetime_str
        
        self.sampling_rate = int(mne_raw.info["sfreq"])
        self.nyquist_frequency = self.sampling_rate/2
        self.channels = mne_raw.info["ch_names"]

        self.online_high_pass = self.info["highpass"]
        self.online_low_pass = self.info["lowpass"]
        if self.online_high_pass == 0:
            self.online_high_pass = None #This corresponds to a DC signal being recorded (i.e., no online filter being applied)
        if self.online_low_pass == self.sampling_rate/2:
            self.online_low_pass = None #Per MNE raw documentation, if no low pass filter is applied online, then self.info["lowpass"] will record half the sampling rate.
        
        channels_dict = {'CrossEeg': 'EEG Cross'} #Continue adding as needed
        if any(channel not in self.channels for channel in which_channel):
            for i in range(len(which_channel)):
                    which_channel[i] = channels_dict[which_channel[i]]


        if any(channel not in self.channels for channel in which_channel):
            for i, channel in enumerate(self.channels):
                print(f"{i + 1}: {channel}")
            raise ValueError(f"You entered a channel name in the which_channel variable that is not in the MNE raw object. The specific channel(s) that you wanted was/were: {[channel for channel in which_channel if channel not in self.channels]}. The channels available are {self.channels}.")
        
        self.channels=which_channel
        self.viewTime1 = viewTime1
        self.viewTime2 = viewTime2

        if len(which_channel) == 1:
            if which_channel[0] not in self.channels: #Forces user to choose from list of channels if the input one does not exist
                which_channel = [choose_channel(self)]
                self.channels = which_channel
            self.EEG_data, self.EEG_times_from_zero = mne_raw.get_data(picks = which_channel,return_times = True, units='uV') #VERY IMPORTANT; self.EEGtimes_from_zero are in SECONDS
            if isinstance(self.EEG_data, np.ndarray):
                # If it is, convert to a regular Python list
                EEGdata = self.EEG_data
                if len(EEGdata)==1:
                    self.EEG_data = EEGdata.tolist()[0]
                else:
                    pass
            if isinstance(self.EEG_times_from_zero, np.ndarray):
                # If it is, convert to a regular Python list
                EEGtimes = self.EEG_times_from_zero
                if len(EEGtimes)==1:
                    self.EEG_times_from_zero = EEGtimes.tolist()[0]
                else:
                    pass
        elif len(which_channel) == 2: #This is how we will build bipolar signals
            mne_raw.get_data(picks = which_channel,return_times = True, units='uV')
            EEG_data_arr, self.EEG_times_from_zero = mne_raw.get_data(picks = which_channel,return_times = True, units='uV')
            self.EEG_data = EEG_data_arr[0]-EEG_data_arr[1]
            
        def __eq__(self, other): 
            return self.__dict__ == other.__dict__

    @property
    def signal_length_secs(self):
        """Recalculate and return the signal length in seconds."""
        if len(self.EEG_times_from_zero) == 0:
            return 0  # Handle the case where EEG_times_from_zero might be empty
        return self.EEG_times_from_zero[-1]  # Return the last value in the array

    @property
    def chunks(self):
        return self.windows

    @property
    def file_name_whole(self):
        my_str = os.path.splitext(os.path.basename(self.file_path_current))[0] #will return the entire file name, just without .pkl (e.g., FIS101post_P3-P4_stitched_cut_artefact_rejected_20_24_hrs)
        return my_str

    @property
    def file_name(self):
        my_str = self.file_name_whole.split('_')[0] #Will return FIS102, for example
        return my_str
    
    @property
    def fis(self):
        base_str = self.file_path_original
        # Look for FIS followed by digits between backslashes
        match = re.search(r'\\(FIS\d+)', base_str, re.IGNORECASE)
        if match:
            fis_part = match.group(1).upper()
            return fis_part
        else:
            raise ValueError("This file's FIS could not be extracted from its file name.")
        
    @property
    def num_fis(self):
        base_str = self.fis
        return ''.join(char for char in base_str if char.isdigit()) #(i.e., 102 from a FIS102 file)
    
    @property
    def analyze_end_time(self):
        return self.rec_start_date_time
    
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

    def get_val():
        """
        Function with the goal of finding a value in a given Excel file
        """

    def mne_neuralPy_compare(self,bool_visualize = False):
        
        """
        This method simply checks that the EEG voltages extracted by the neuralPy methods and MNE methods are identical. NB: we are forcing both to examine C3, since this is just an exercise.

        \nTakes in:
        bool_visualize (boolean): If True, we will get plots of EEG signals from the neuralPy and MNE methods.

        \nReturns:
        boolean_identical: A boolean value that is true if the EEG voltage arrays are identical, false if they are not
        """
        
        neuralPy_data_out, ch_fs = preprocessing_EEG.load_edf_EEG_file(self.file_path,channels_to_look=['C3','C4','P3', 'P4'])
        neuralPy_data = neuralPy_data_out.loc[self._s1:self._s2, 'C3']
        MNE_data = self.EEG_data

        if bool_visualize:
            #Plot neuralPy data
            plt.figure(figsize=(20, 5))
            plt.plot(neuralPy_data, 'b-', label='NeuralPy_C3') #This plots x,y as time values, voltage values
            plt.grid(linestyle='-', linewidth=1)
            plt.legend()
            plt.show()

            #Plot MNE data
            plt.figure(figsize=(20, 5))
            plt.plot(self.EEG_times_from_zero, MNE_data, 'r-', label='MNE_C3')
            plt.grid(linestyle='-', linewidth=1)
            plt.legend()
            plt.show()

        count = 0
        countWrong = 0
        for i, j in zip(MNE_data, neuralPy_data):
            if i != j:
                print(str(i) + " " + str(j))
                countWrong += 1
            count += 1
        print("Across "+str(count)+" items, "+str(countWrong)+" did not match.")
        if countWrong == 0:
            return True
        else:
            return False
        #Because this returns true, we know that EEG C3-REF in the MNE raw object is the same as what we get from data_out on NeuralPy. Since the documentation in MNE is more robust, I will continue there  """
    def add_change_log(self, changeText, suppress=False):
        """
        This method is used to write to an object's change-log every time something is done. This was unfortunately not implemented into the functions themselves, so this method must be called AFTER EACH MODIFICATION WE WANT TO LOG.

        Inputs:
        - changeText: the text (as a string) that will be added into the changelog

        Outputs: none (this method only updates our changelog)
        """
        now = datetime.now()
        datetime_str = now.strftime("%Y-%m-%d %H:%M:%S")
        if suppress == False:
            print(datetime_str+"\n"+changeText)

        self.change_log = self.change_log + "\n\n"+datetime_str+"\n"+changeText

    def save_as(self, save_directory, suffix):
        """
        This method is used to save an alex_raw object to a given directory. This will be useful once a full processing pipeline is built, so that the object can be checked up on at each interpretation stage.

        Inputs:
        - save_directory (str): the path directory where we would like the file to be saved
        - suffix (str): any text we would like appended to the end of the file name (used for keeping track of which processing step each file was used for). NB: we DO NOT have a way to edit the path currently

        Outputs: none (this method just saves the file)
        """
        if ".pkl" not in suffix:
            suffix = suffix+".pkl"
        else:
            pass
        
        #Update the file path with the suffix we want; self.file_path has the file path we're interested in
        if self.saved_once == False: #The file name updating is done only the first save
            path = self.file_path_current
            directory, filename = os.path.split(path) #For dir/subdir/filename.extension, returns dir/subdir/ as directory and filename.extension as extension
            directory = save_directory #This allows us to build a new directory just for saved files
            name, extension = os.path.splitext(filename) #for filename.extension, returns filename as name and .extension as extension
            name = name + "_" + self.channel_name.replace(" ", "_")
            #extension = ".pkl" #When we save our files using the dill module, we no longer have a readable .edf file, but instead something new
            filename_new = name+suffix #This is the new file name of the object we will save
            filepath_current = os.path.join(directory, filename_new)
            self.file_path_current = filepath_current #Our object's current path is updated (as this is the variable that will be used in future iterations of saveAs). The original file lxn/name can be accessed in the variable self.file_path_original
            self.file_dir_current = directory
            self.saved_once = True

        self.add_change_log(f"File saved to {self.file_path_current}")
        with open(self.file_path_current, 'wb') as file:
            pickle.dump(self, file, protocol=pickle.HIGHEST_PROTOCOL)
        
        import gc
        gc.collect()
        
        print(f"File {self.file_path_current} saved.")

    def save(self):
        """
        This method is used to save an alex_raw object that has already been saved once (i.e., that already has a directory and file_name established) as its same file name to the same file_directory.

        Inputs:
        None (all information we need is contained in the self object).

        Outputs: none (this method just saves the file)
        """
        if self.saved_once == False:
            def select_directory():
                # Initialize the Tkinter root widget
                root = tk.Tk()
                root.withdraw()  # Hide the root window

                # Open a dialog to select a directory
                root.attributes('-topmost', True)  # Bring the dialog to the front
                directory_path = filedialog.askdirectory(title="Select a directory")
                
                # Destroy the root window after selection
                root.destroy()
                
                # Return the selected directory path
                return directory_path
            directory_path = select_directory()
            self.save_as(directory_path, ".pkl")
        else:
            self.add_change_log(f"File saved to {self.file_path_current}")
            with open(self.file_path_current, 'wb') as file:
                pickle.dump(self, file, protocol=4)

    def downsample_EEG(self, times, data, new_sample_frequency):
        """
        Downsample EEG data to a lower sampling rate.
        
        Parameters:
        EEG_times (array): Array of time points corresponding to the EEG data.
        EEG_data (array): Array of EEG voltage data.
        original_sampling_rate (float): Original sampling rate of the data (Hz).
        downsampling_rate (float): Desired sampling rate (Hz).
        
        Returns:
        downsampled_times (numpy array): Downsampled time points.
        downsampled_data (numpy array): Downsampled EEG data.
        """
       
        if new_sample_frequency == None:
            return times, data

        downsampled_times = []
        downsampled_data = []
        original_srate = int(self.sampling_rate) #Eg 250 (per sec)
        interval = int(original_srate//new_sample_frequency) #If we want 10 samples per sec, this becomes 250/10 = every 25th data point. Always rounds down so we never go over a given second in sampling each second. We built the loop to reset every second (and thereby avoid drift with the rounding over very large data sets)
        downsampled_times = times[::interval]
        downsampled_data = data[::interval]
        reduction_factor = len(data)/len(downsampled_data)
        print(f"Downsampled data. Data reduced by a factor of {reduction_factor}")
        return downsampled_times, downsampled_data
    
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
            windows_hr = [[(elem / 3600.0) for elem in win] if isinstance(win, (list, tuple)) else win / 3600.0
                for win in windows_in_seconds]
            start_time_hr = round(self.EEG_times_from_analyze_end_time[0] / 3600, 2)
            end_of_surgery = self.analyze_end_date_time.strftime("%d %B, %Y %H:%M")
            signal_duration_hr = round(self.signal_length_secs / 3600, 2)
            raise WindowBuildError(f"None of this file was in the windows of interest (hr): {windows_hr}. The file started {start_time_hr} hours after the end of surgery (which was on {end_of_surgery}) and ran {signal_duration_hr} hrs.")
        else:
            return indices_arr, windows_arr

    def view(self, view_time_1=None, view_time_2=None, view_ind_1 = None, view_ind_2 = None, downsample = None, view_data=None, view_times=None, side=None, block_bool=True, plt_title = None):
        """
        Function that displays EEG signal at the time we're interested in.

        Inputs:
        - view_viewTime1 (int): time in the view_times array we want to start viewing
        - view_viewTime2 (int): time in the view_times array we want to stop viewing
        - downSample (int): the frequency we would like to downsample to (currently shows at sfreq)
        - view_data (arr): the data we will be viewing
        - view_times (arr): the time points for the data we will be viewing
        - side (str): if not None, can be "left" or "right"; this will put the output figure on the left or right side
        - block_bool (bool): sets block variable in our plot. If True, code will not proceed until the plot is closed.
        - plt_title (str): title of the plot
        """
        print("Viewing data.")
        if not view_times is None:
            if len(view_times)==0:
                print("No view times given, so this data cannot be viewed.")
                return 0
        else:
            view_times = self.EEG_times_from_zero
        if view_data is None:
            view_data = self.EEG_data

        view_times = np.array(view_times)
        view_data = np.array(view_data)
        if not(len(view_times)==len(view_data)):
            raise ValueError(f"Your view_time variable (length {len(view_times)} samples) is not the same length as the data array we're using (length {len(view_data)})")

        if downsample is None:
            downsample = self.sampling_rate
        if view_time_1 is None:
            view_time_1 = self.viewTime1
        elif view_time_1<0:
            view_time_1 = self.viewTime1
            print("Initial view time was under zero, so corrected to zero")
        if view_time_2 is None:
            view_time_2 = self.viewTime2
        elif view_time_2>view_times[-1]:
            print("Final view time was larger than signal time. Will show until end of signal.")
            view_time_2 = view_times[-1]
        if view_data is None:
            view_data = self.EEG_data
        if view_times is None:
            view_times = self.EEG_times_from_zero
            
        print("Determining view indices from the data")
        view_index_1 = np.where(view_times >= view_time_1)[0][0]
        view_index_2 = np.where(view_times <= view_time_2)[0][-1]
        if not view_ind_1 is None:
            view_index_1 = view_ind_1
        if not view_ind_2 is None:
            view_index_2 = view_ind_2

        view_times = view_times[view_index_1:view_index_2]
        view_data = view_data[view_index_1:view_index_2]
        #view_times, view_data = self.downsample_EEG(view_times, view_data, downsample)

        #Turns NaN values into zeros JUST FOR VIEWING
        if np.isnan(view_times).any():
            view_data = np.nan_to_num(view_data, nan=0)


        def move_existing_plot_to_half_screen(side="left"):
            
            """
            Inputs:
            - side (str): can be "left" or "right"; describes what side of the screen the plot should show on
            """
            
            # Find the plot window by title (usually 'Figure' for matplotlib plots)
            try:
                plot_window = gw.getWindowsWithTitle('Figure')[0]
            except IndexError:
                print("No plot window found. Make sure the plot is open.")
                return

            # Get screen dimensions
            screen_width, screen_height = pyautogui.size()

            # Calculate new window size (half the screen width, full screen height)
            window_width = screen_width // 2
            window_height = screen_height

            # Set window position and size
            if side=="left":
                # Move to the left half of the screen
                plot_window.moveTo(0, 0)
            elif side=="right":
                # Move to the right half of the screen
                plot_window.moveTo(window_width, 0)

            # Resize the window to take up half the screen width and full height
            plot_window.resizeTo(window_width, window_height)


        plt.figure(figsize=(20, 5))
        file_name = os.path.basename(self.file_path_current)
        plt.plot(view_times, view_data, color='blue', label=f'EEG signal at channel: {self.channel_name} in the file {file_name}')
        # Add grid, legend, and show plot
        plt.grid(linestyle='-', linewidth=0.25)
        #plt.title("Original signal seen from " + str(viewTime1/250) + " sec to " + str(viewTime2/250) + "sec")
        plt.legend(loc="lower left")
        if not side==None:
            move_existing_plot_to_half_screen(side)
        else:
            manager = plt.get_current_fig_manager()
            manager.window.geometry("1366x768-30+30")
        if not plt_title is None:
            plt.title(plt_title)
        else:
            plt.title("Current EEG data")
        plt.show(block=block_bool) #allows code to continue to run without closing plot

    def view_filtered(self, filtered_data, original_data=None, view_viewTime1=None, view_viewTime2=None, downsample = None, view_times=None, block_bool=False, plt_title=None, plot_power_spectrum = True, power_spectrum_title = None):
        """
        This is a function meant to be used when an alex_raw object is being filtered and we would like to see the the filtered signal superimposed on the raw one.

        Inputs:
        - view_viewTime1 (int or float): time in the view_times array we want to start viewing
        - view_viewTime2 (int or float): time in the view_times array we want to stop viewing
        - downSample (int): the frequency we would like to downsample to (default is sfreq)
        - view_data (arr): the data we will be viewing
        - view_times (arr): the time points for the data we will be viewing
        - block_bool (boolean): if True, will prevent code from advancing until the displayed plot is closed
        - plt_title (str): title of the plot

        Outputs:
        - None (simply outputs the chart)
        """
        #Initialize variables as needed
        downsamp_bool = True
        if downsample is None:
            downsample = self.sampling_rate
            downsamp_bool = False
        if view_viewTime1 is None:
            view_viewTime1 = self.viewTime1
        elif view_viewTime1<0:
            view_viewTime1 = self.viewTime1
            print("Initial view time was under zero, so corrected to zero")
        if view_viewTime2 is None:
            view_viewTime2 = self.viewTime2
        elif view_viewTime2>view_times[-1]:
            print("Final view time was larger than signal time. Will show until end of signal.")
            view_viewTime2 = view_times[-1]
        if original_data is None:
            original_data = self.EEG_data
        if view_times is None:
            view_times = self.EEG_times_from_zero
        if not (len(original_data) == len(view_times) and len(filtered_data)==len(view_times)):
            raise ValueError("The arrays of EEG data and EEG time points are not the same length. Please review your inputs for view_data and view_time.")
        print("Determining view indices from the data")
        # Find the index closest to view_viewTime1
        import bisect
        view_index_1 = bisect.bisect_left(view_times, view_viewTime1)
        if view_index_1 == len(view_times):  # If it's beyond the last element
            view_index_1 -= 1
        elif view_index_1 > 0 and abs(view_times[view_index_1] - view_viewTime1) >= abs(view_times[view_index_1 - 1] - view_viewTime1):
            view_index_1 -= 1
        # Find the index closest to view_viewTime2
        view_index_2 = bisect.bisect_left(view_times, view_viewTime2)
        if view_index_2 == len(view_times):  # If it's beyond the last element
            view_index_2 -= 1
        elif view_index_2 > 0 and abs(view_times[view_index_2] - view_viewTime2) >= abs(view_times[view_index_2 - 1] - view_viewTime2):
            view_index_2 -= 1
        view_times = view_times[view_index_1:view_index_2]
        original_data = original_data[view_index_1:view_index_2]
        if np.isnan(original_data).any():
            original_data = np.nan_to_num(original_data, nan=0)
        filtered_data = filtered_data[view_index_1:view_index_2]
        if np.isnan(filtered_data).any():
            filtered_data = np.nan_to_num(filtered_data, nan=0)
        if downsamp_bool:
            view_times, original_data = self.downsample_EEG(view_times, original_data, downsample)
            view_times, filtered_data = self.downsample_EEG(view_times, filtered_data, downsample)
        
        #Visualize
        plt.figure(figsize=(20, 5))
        plt.plot(view_times[view_index_1:view_index_2], original_data[view_index_1:view_index_2], color='blue', label='Original EEG signal at channel: '+self.channel_name)
        # Plot the filtered signal in red
        plt.plot(view_times[view_index_1:view_index_2],filtered_data[view_index_1:view_index_2], color='#FF8C00', linewidth=2, label='Filtered EEG signal at channel: '+self.channel_name)
        # Add grid, legend, and show plot
        plt.grid(linestyle='-', linewidth=1)
        #plt.title("Original signal seen from " + str(viewTime1/250) + " sec to " + str(viewTime2/250) + "sec")
        plt.legend()
        if not plt_title is None:
            plt.title(plt_title)
        else:
            plt.title("Original and Filtered EEG data")
        plt.show(block = block_bool)

        if plot_power_spectrum == True:
            # Plot raw data power spectrum (in black)
            frequencies_raw, psd_raw = welch(original_data[:], fs=self.sampling_rate, nperseg=1024)
            # Plot new (filtered) power spectrum (in green)
            frequencies_filt, psd_filt = welch(filtered_data[:], fs=self.sampling_rate, nperseg=1024)  # Use only the first EEG channel (C3-REF)
            # Create a plot to display both PSDs
            plt.figure(figsize=(10, 6))
            # Plot raw data power spectrum in black
            plt.semilogy(frequencies_raw, psd_raw, 'k-', label='Raw EEG at channel ' +self.channel_name)  # 'k-' for black line
            # Plot filtered data power spectrum in green
            plt.semilogy(frequencies_filt, psd_filt, color="orange", label='Filtered EEG at channel '+self.channel_name)  # 'g-' for green line
            # Add plot details
            if not power_spectrum_title is None:
                plt.title(power_spectrum_title)
            else:
                plt.title('Power Spectrum of Raw and Filtered EEG Signal'+self.channel_name)
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('Power Spectral Density (dB/Hz)')
            plt.xlim([0, 60])  # Plot up to 60 Hz (adjust as needed)
            plt.grid(True)
            plt.legend()  # Show legend to differentiate between raw and filtered
            # Display the plot
            plt.show()


    def view_rEEG(self, rEEG_data, rEEG_times, original_data=None, original_times=None, view_time_1=None, view_time_2=None, view_original = True, block_bool=False, plt_title=None, n = 1):
        """
        This is a function meant to be used when an alex_raw object is being filtered and we would like to see the the filtered signal superimposed on the raw one.

        Inputs:
        - rEEG_data (arr of float): rEEG data of interest
        - rEEG_times (arr of float): time points IN SECONDS corresponding to our rEEG data
        - original_data (arr): the original data we will be viewing
        - original_times (arr): time points IN SECONDS from the original data
        - view_time_1 (int or float): time point IN SECONDS where we want to start viewing
        - view_time_2 (int or float): time point IN SECONDS where we want to stop viewing
        - view_original (boolean): if True, will display the original signal with the rEEG signal
        - block_bool (boolean): if True, will prevent code from advancing until the displayed plot is closed
        - plt_title (str): title of the plot
        - n (int or float): our viewing fudge factor; we divide view_time_2 by this factor to change the resolution of the viewing plot. If n is 3, for example, view_time_2 becomes view_time_1 + 1/3(the original view time); that is, we only see the first third of the signal

        Outputs:
        - None (simply shows the plot)
        """
        #Initialize variables as needed
        if original_data is None:
            original_data = self.EEG_data
        original_data = np.where(np.isnan(original_data), 0, original_data) #Replace NaNs with 0s for viewing
        
        if original_times is None:
            original_times = self.EEG_times_from_analyze_end_time
        if not len(original_data) == len(original_times):
            raise ValueError("The arrays of original EEG data and original EEG time points are not the same length. Please review your inputs for original_data and original_times.")
        
        #CONTINUE CODING HERE
        if not (len(rEEG_data) == len(rEEG_times)):
            raise ValueError("The arrays of EEG data and EEG time points are not the same length. Please review your inputs for view_data and view_time.")
        if view_time_1 is None:
            view_time_1 = rEEG_times[0]
        elif view_time_1<0:
            view_time_1 = rEEG_times[0]
            print("Initial view time was under zero, so corrected to zero")
        if view_time_2 is None:
            view_time_2 = rEEG_times[-1]
        elif view_time_2>rEEG_times[-1]:
            print("Final view time was larger than the original signal time. Will show until end of original signal.")
            view_time_2 = rEEG_times[-1]
        if not n==1:
            view_time_2 = round(view_time_1+(view_time_2-view_time_1)/n,1) #

        print("Determining view indices from the data")
        # Find the index of the value closest to (but not under) view_time_1 in rEEG_times
        rEEG_index_start = np.searchsorted(rEEG_times, view_time_1, side="left")
        rEEG_index_end = np.searchsorted(rEEG_times, view_time_2, side="right")
        original_index_start = np.searchsorted(original_times, view_time_1, side="left")
        original_index_end = np.searchsorted(original_times, view_time_2, side="right")

        
        #Visualize
        if view_original:
            plt.figure(figsize=(20, 5))
            
            plt.plot(original_times[original_index_start:original_index_end], original_data[original_index_start:original_index_end], color='blue', label=f'Original EEG signal for file {self.file_name} at channel: {self.channel_name}')
            plt.show(block=False)

        plt.plot(rEEG_times[rEEG_index_start:rEEG_index_end],rEEG_data[rEEG_index_start:rEEG_index_end], color='#FF8C00', linewidth=2, label=f'rEEG signal for file {self.file_name} at channel: {self.channel_name}')
        plt.grid(linestyle='-', linewidth=1)
        plt.legend()
        if not plt_title is None:
            plt.title(plt_title)
        else:
            plt.title(f"rEEG data for file {self.file_name}")
        plt.show(block = block_bool)

    
    def save_as_easy(self, directory, suffix):
        """
        Function that helps us save a file while adding suffix. self.save_as has been optimized to work with the EDF reading step, whereas this function is good for everything else.

        Inputs:
        - directory (str) = the path where we would like the file to be dropped 
        - suffix (str): the suffix we would like to add onto the end of the file name (but before the .pkl file ending)
        """
        current_path = self.file_path_current
        file_name_no_pkl = os.path.splitext(os.path.basename(current_path))[0]
        new_file_name = file_name_no_pkl+suffix+'.pkl'
        new_path = os.path.join(directory, new_file_name)
        self.file_path_current = new_path
        self.add_change_log(f"File saved to {new_path}",suppress=True)
        with open(new_path, 'wb') as file:
            pickle.dump(self, file, protocol=pickle.HIGHEST_PROTOCOL)
        import gc
        gc.collect()
        print(f"File {self.file_path_current} saved.")