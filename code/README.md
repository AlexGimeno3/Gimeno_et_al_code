# General info:
This repository is associated with the paper [paper title, etc here].

# About analyze_files.py
***Please note: this repository is intended to allow for reviewers/other researchers to examine the methodology in the paper [TITLE]. It has not been validated for clinical use. Please DO NOT use it for any medical purposes.***

This code was designed to import raw EEG files (specifically, .edf files), process them, and then build a GUI to allow for tailored analysis. The purpose of this document is to provide the reader with more methodological detail than is possible in the main manuscript or supplementary materials. Details for each .py file are provided at the top of each file.

# Pre-analysis
## Downloading EDF files
EEG files were downloaded as described in the PDFs in the EEG_file_downloading directory. The options _Entire signal_ and _Raw data_ were selected during download. The file readers on the Nicolet machines included in the final sample were:  
- NicoletOne Monitor Reader 5.71.0.2439
- NicoletOne Monitor Reader 5.94.1.534

## EDF file storage
This code was designed based on the specific file formats and structures that were encountered using Nicolet and OBM (Olympic Brainz Monitor) machines; therefore, the expected directory structures of the EEG files are built around these machines. (NB: OBM machine data was omitted from final analysis, but as they are included in the current code, their storage/processing has been maintained for completeness. The code is built such that OBM files can be analyzed as well, should the user select this.)

### Nicolet
All Nicolet folders were stored in one subdirectory. A given FIS (that is, ID) may have multiple folders; however, each folder contained the EDF files from one single recording. Folder names began with "FISXXX_" (where XXX is any integer over 0). Within the folder, each edf file started with "FISXXX". If a single recording was split amongst multiple EDF files (this happens when EDF-D files are encountered/exported on a machine), the first had no suffix, followed by 1, then 2, then 3, etc.

### OBM
All OBM folders were stored in one subdirectory. A given FIS may have had multiple folders (each corresponding to a different recording on the same machine). If just one folder for a given FIS was present, the folder's name was "FISXXX"; if there were multiple, the folders' names were "FISXXXa", "FISXXXb", "FISXXXc", etc. The current code was built to handle the case where, in a given FIS's folders, there were numerous subdirectories, the key ones being called "CrossEeg", "LeftEeg" and "RightEeg" (note how FIS number was not mentioned in these folders). Within each was found at least 2 EDF files: an original (which contained simply the folder name + .edf; e.g., "CrossEeg.edf") and at least one exported EDF file. The exported EDF files had the name of f"{fis_folder_name}_{subfolder_name}_{n}", where n is always an integer, starting at 1, expressed with 4 digits. For example, if we are looking at the first exported EDF in the CrossEeg folder of the folder FIS92b, the file name of that EDF file must be "FIS92b_CrossEeg_0001.edf". The next would be "FIS92b_CrossEeg_0002.edf", etc. Exportation from the EDF-D file (i.e., "CrossEeg.edf") was done using EDFBrowser version 2.12 64-bit, specifically the "Convert EDF+D to EDF+C" tool.

## Other necessary files:
A few more data files were needed to aid with determining recording times relative to surgery start and end times.
### Nicolet recording times
This Excel file was called "nicolet_file_times.xlsx" and contained the following columns: "FIS" (with entries as strings such as FIS1), "fis_num" (with entries as integers such as 1), "folder_path" (the path to the folder containing the downloaded FIS file), "date_cx" (the date of the surgeries), "analyze_end_time" (the time that the surgery ended, expressed as HHMM), "recording_start_date" (the date the recording started, expressed as DDMMYYY), "recording_start_time" (expressed as HHMMSS), and "recording_time_total" (expressed as HHMMSS).
### OBM recording times
This Excel file was called and "OBM_file_times.xlsx" and contained the following columns: "folder_path" (the path to the folder containing the downloaded FIS file), "recording_start_date" (the date the recording started, expressed as DDMMYYY), and "recording_start_time" (expressed as HHMMSS).
### Surgery times
**NB: This code has neither been built for nor validated on surgeries ocurring during two different calendar days (i.e., starting at 2300 on Jan 1 and ending at 0100 on Jan 2).**
This Excel file was named "surgery_time_data.xlsx" and contained the following columns: "fis_num" (with entries as integers such as 1), "date_surgery" (with entries as dates formatted as 1900-04-16), "time_start_cx" (with entries as surgery start times formatted as HHMM), "time_end_cx" (with entries as surgery end times formatted as HHMM), "time_start_ecc" (with entries as the start of CPB formatted as HHMM), "time_start_clamp" (with entries as the start of aortic cross clamping formatted as HHMM), "time_end_clamp" (with entries as the end of aortic cross-clamping formatted as HHMM), and "time_end_ecc" (with entries as the end of CPB formatted as HHMM).

## Dependencies
The analysis was run on Python 3.12.7, 64-bit. 
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
