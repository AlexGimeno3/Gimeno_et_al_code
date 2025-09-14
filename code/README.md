# General info
This directory contains all the files used to complete the analysis in the paper _Association of quantitative EEG parameters extracted during infant cardiopulmonary bypass and two-year neurodevelopmental outcomes_.

Each subdirectory contains a README file as necessary to explain the code and structure of the pipeline.

# About analyze_files.py
***Please note: this repository is intended to allow for reviewers/other researchers to examine the methodology used in the paper _Association of quantitative EEG parameters extracted during infant cardiopulmonary bypass and two-year neurodevelopmental outcomes_. It has not been validated for clinical use. Please DO NOT use it for any medical purposes.***

This is the "master script" that, when run, will perform all the the analyses needed to duplicate the analysis performed in this paper. Note that it contains one setting to change in the "__main__" section at the bottom of the script: "postop" or "intraop". This setting directs the script to either access "processing_parameters_intraop.py" or "processing_parameters_postop.py" (in the ./utils subdirectory). **NB: to edit the analysis parameters, these processing_parameters files are the ones to alter.**

# Dependencies
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
