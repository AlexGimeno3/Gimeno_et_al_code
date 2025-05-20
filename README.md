# General info:
This repository is associated with the paper [paper title, etc here].

It is organized as follows: [organization here]


# About run_all.py
***Please note: this program is intended to allow for reviewers/other researchers to examine the methodology in the paper [TITLE]. It has not been validated for clinical use. Please DO NOT use it for such.***

This code is designed to take in raw EEG files (specifically, .edf files), process them, and prepare all necessary objects for application of stats to the sample. The purpose of this document is to walk the user through the processing steps in each script as well as the necessary dependencies. Please read this in its entirety before using the code.

# Before starting
## Downloading EDF files
NB: While downloading, it would be a smart idea to note recordin start and end dates and times (see the section **Other Necessary Files**)
### For Nicolet machines
[Add Nicolet details here, including how to download.]

### For OBM machines
[Add OBM details here, including how to download].

After this downloading procedure, I recommend, for each FIS folder, the user programmatically create a subdirectory for each downloaded EDF file (exactly sharing that file's name, except the suffix) and move that EDF file into its respective new subdirectory. This will speed up the next step in OBM file preparation (see next paragraph).

The OBM files will be downloaded in EDF-D format (i.e., as discontinuous files); MNE cannot read EDF-D files, only EDF-C files. Therefore, the user will need to manually convert EDF-D to EDF-C files. I opted to do this using EDFbrowser (https://www.teuniz.net/edfbrowser/) version 2.12, 64-bit. Specifically, the user can go to Tools --> Convert EDF-D to EDF-C, which will create the necessary continuous EDF files for processing. These files will also contain the necessary metadata to merge into one data stream for later processing.

## Storing EDF files
This code was designed based on the specific file formats and structures that I encountered using Nicolet and OBM machines [GIVE HARDWARE/SOFTWARE details of the EEG machines used]. Therefore, the expected directory structures of the EEG files are built around these machines.

### Nicolet
All Nicolet folders should be stored in one subdirectory. A given FIS (that is, ID) may have multiple folders; however, each folder must have the EDF files from one single recording. Folder names must begin with "FISXXX_" (where XXX is any integer over 0). Within the folder, each edf file must start with "FISXXX". If a single recording is split amongst multiple EDF files (this happens when EDF-D files are encountered/exported on a machine), the first should have no suffix, followed by 1, then 2, then 3, etc.

### OBM
All OBM folders should be stored in one subdirectory. A given FIS may have multiple folders. If just one folder for a given FIS is present, the folder's name should be "FISXXX"; if there are multiple, the folders' names should be "FISXXXa", "FISXXXb", "FISXXXc", etc. The current code is built to handle the case where, in a given FIS's folders, there are numerous subdirectories, the key ones being called "CrossEeg", "LeftEeg" and "RightEeg" (note how FIS number is not mentioned in these folders). Within each should be found at least 2 EDF files: an original (which contains simply the folder name + .edf; e.g., "CrossEeg.edf") and at least one exported EDF file. The exported EDF files should have a name of f"{fis_folder_name}_{subfolder_name}_{n}", where n is always an integer, starting at 1, expressed with 4 digits. For example, if we are looking at the first exported EDF in the CrossEeg folder of the folder FIS92b, the file name of that EDF file must be "FIS92b_CrossEeg_0001.edf". The next would be "FIS92b_CrossEeg_0002.edf", etc. 

## Other necessary files:
A few more data files are needed to aid with determining recording times relative to surgery start and end times.
### Nicolet recording times
This Excel file should contain the following columns: "FIS" (with entries as strings such as FIS1), fis_num (with entries as integers such as 1), folder_path (the path to the folder containing the downloaded FIS file), date_cx (the date of the surgeries), analyze_end_time (the time that the surgery ended, expressed as HHMM), recording_start_date (the date the recording started, expressed as DDMMYYY), recording_start_time (expressed as HHMMSS), and recording_time_total (expressed as HHMMSS).
### OBM recording times
Pass
### Surgery times
_NB: The code has neither been built for nor validated on surgeries ocurring over two days (i.e., starting at 2300 on Jan 1 and ending at 0100 on Jan 2)._ It is very much recommended to omit these surgeries.

Dependencies:
- os v
- sys v
- 
