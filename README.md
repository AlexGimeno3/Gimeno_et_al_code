# General info:
This repository is associated with the paper [paper title, etc here].

# About analyze_files.py
***Please note: this repository is intended to allow for reviewers/other researchers to examine the methodology in the paper [TITLE]. It has not been validated for clinical use. Please DO NOT use it for any medical purposes.***

This code was designed to import raw EEG files (specifically, .edf files), process them, and then build a GUI to allow for tailored analysis. The purpose of this document is to provide the reader with more methodological detail than is possible in the main manuscript or supplementary materials. Details for each .py file are provided at the top of each file.

# Pre-analysis
## Downloading EDF files
### For Nicolet machines
[Add Nicolet details here, including how to download.]

### For OBM machines
[Add OBM details here, including how to download].

After this downloading procedure, I recommend, for each FIS folder, the user programmatically create a subdirectory for each downloaded EDF file (exactly sharing that file's name, except the suffix) and move that EDF file into its respective new subdirectory. This will speed up the next step in OBM file preparation (see next paragraph).

The OBM files will be downloaded in EDF-D format (i.e., as discontinuous files); MNE cannot read EDF-D files, only EDF-C files. Therefore, the user will need to manually convert EDF-D to EDF-C files. I opted to do this using EDFbrowser (https://www.teuniz.net/edfbrowser/) version 2.12, 64-bit. Specifically, the user can go to Tools --> Convert EDF-D to EDF-C, which will create the necessary continuous EDF files for processing. These files will also contain the necessary metadata to merge into one data stream for later processing.

## Storing EDF files
This code was designed based on the specific file formats and structures that I encountered using Nicolet and OBM (Olympic Brainz Monitor) machines. Therefore, the expected directory structures of the EEG files are built around these machines. Specific details of the machines are:

- OBM viewer 3.1.5.57
- NicoletOne Monitor Reader 5.71.0.2439
- NicoletOne Monitor Reader 5.94.1.534

### Nicolet
All Nicolet folders should be stored in one subdirectory. A given FIS (that is, ID) may have multiple folders; however, each folder must have the EDF files from one single recording. Folder names must begin with "FISXXX_" (where XXX is any integer over 0). Within the folder, each edf file must start with "FISXXX". If a single recording is split amongst multiple EDF files (this happens when EDF-D files are encountered/exported on a machine), the first should have no suffix, followed by 1, then 2, then 3, etc.

### OBM
All OBM folders should be stored in one subdirectory. A given FIS may have multiple folders. If just one folder for a given FIS is present, the folder's name should be "FISXXX"; if there are multiple, the folders' names should be "FISXXXa", "FISXXXb", "FISXXXc", etc. The current code is built to handle the case where, in a given FIS's folders, there are numerous subdirectories, the key ones being called "CrossEeg", "LeftEeg" and "RightEeg" (note how FIS number is not mentioned in these folders). Within each should be found at least 2 EDF files: an original (which contains simply the folder name + .edf; e.g., "CrossEeg.edf") and at least one exported EDF file. The exported EDF files should have a name of f"{fis_folder_name}_{subfolder_name}_{n}", where n is always an integer, starting at 1, expressed with 4 digits. For example, if we are looking at the first exported EDF in the CrossEeg folder of the folder FIS92b, the file name of that EDF file must be "FIS92b_CrossEeg_0001.edf". The next would be "FIS92b_CrossEeg_0002.edf", etc. 

## Other necessary files:
A few more data files are needed to aid with determining recording times relative to surgery start and end times.
### Nicolet recording times
This Excel file should contain the following columns: "FIS" (with entries as strings such as FIS1), "fis_num" (with entries as integers such as 1), "folder_path" (the path to the folder containing the downloaded FIS file), "date_cx" (the date of the surgeries), "analyze_end_time" (the time that the surgery ended, expressed as HHMM), "recording_start_date" (the date the recording started, expressed as DDMMYYY), "recording_start_time" (expressed as HHMMSS), and "recording_time_total" (expressed as HHMMSS).
### OBM recording times
This Excel file should contain the following columns: "folder_path" (the path to the folder containing the downloaded FIS file), "recording_start_date" (the date the recording started, expressed as DDMMYYY), and "recording_start_time" (expressed as HHMMSS).
### Surgery times
***NB: The code has neither been built for nor validated on surgeries ocurring during two different calendar days (i.e., starting at 2300 on Jan 1 and ending at 0100 on Jan 2).*** It is very much recommended to omit these surgeries.

This Excel file should contain the following columns: "fis_num" (with entries as integers such as 1), "date_surgery" (with entries as dates formatted as 1900-04-16), "time_start_cx" (with entries as surgery start times formatted as HHMM), "time_end_cx" (with entries as surgery end times formatted as HHMM), "time_start_ecc" (with entries as the start of CPB formatted as HHMM), "time_start_clamp" (with entries as the start of aortic cross clamping formatted as HHMM), "time_end_clamp" (with entries as the end of aortic cross-clamping formatted as HHMM), and "time_end_ecc" (with entries as the end of CPB formatted as HHMM).

## Dependencies
The analysis for this paper was run on Python 3.12.7, 64-bit. 
Dependencies:
 - numpy 1.26.4
 - scipy 1.11.4
 - pandas 2.2.3
 - matplotlib 3.10.0
 - mne 1.9.0
 - neurokit2 0.2.11
 - pyEDFlib 0.1.40
 - pyGetWindow 0.0.9
 - pyAutoGUI 0.9.54
 - openpyxl 3.1.5
 - nolds 0.6.1
 - statsmodels 0.14.4
 - EntropyHub 2.0
 - PyPDF2 3.0.1
 - NEURAL_py_EEG 0.1.4 (NB: I also have a locally downloaded and managed/edited NEURAL_py_EEG folder in utils. This folder should be left untouched, while the dependency here would also need downloaded.)

# Running the code
## Adjusting parameters
In analyze_files.py, the user can designate if they would like to process intraoperative or postoperative data. Before running the script, the user will need to do 2 things:
  1. Update the setting variable on line 83. Valid values are "intraop", if the user is processing intraoperative data, and "postop", if the user is processing postoperative data. Respectively, these point to processing_parameters_intraop.py and processing_parameters_postop.py in the utils folder.
  2. Update processing_parameters_intraop.py or processing_parameters_postop.py. The default values in these dictionaries are the ones that were used for the analysis in [paper title, etc here]. It is recommended that the user not touch the variables under #File structure variables. The other variables can be changed to customize the analysis; descriptions of these variables are included in the respective dictionaries
