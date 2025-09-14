def run_all(time_frame):
    """
    Function to do an entire analysis for us from start to finish.

    Inputs:
    - base_dir (str): the directory in which we will be building all our other files
    - windows (arr): at this point, should either be [["surgery"]] or [[a,b]], where a and b are the start and end times (in HOURS) of the windows of interest.
    """
    import os
    import sys
    # Get the current script's directory (where the current script is located)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.append(os.path.join(current_dir, 'utils'))
    from utils import config
    config.initialize_config(time_frame)
    sys.path.append(os.path.join(current_dir, 'processing_scripts'))
    from _1_nicolet_stitching_cutting import iterate_stitch_n
    from _1_obm_stitching_cutting import iterate_stitch_o
    from _1_all_file_combiner import combine_files
    from _2_artefact_rejection import iterate_artefact_rejection
    from _3_preprocessing import iterate_preprocess
    from _4_extract_variables import extract_vars_iterate
    sys.path.append(os.path.join(current_dir, 'processing_scripts',"5_run_stats"))
    from run_stats import stat_run

    config.initialize_config(time_frame)
    vars_dict = config.get_vars_dict()
    
    windows = vars_dict["windows"]
    if len(windows)>1:
        raise ValueError("The program no longer supports more than one window at a time; however, it is much faster and more memory-efficient. Please make the window a single time range encompassing all the windows of interest. You can then use the GUI to choose the precise window of interest.")
    
    if windows == [["surgery"]]:
        base_dir = os.path.join(vars_dict["base_dir"], "intraop_files")
    else:
        base_dir = os.path.join(vars_dict["base_dir"], "postop_files")
    
    
    #Step 1: Stitch
    print("Stitching...")
    base_dir_step1 = os.path.join(base_dir , "1_stitched_cut")
    iterate_stitch_n(base_dir_step1)
    iterate_stitch_o(base_dir_step1)
    combine_files(base_dir_step1, windows)
    #Metadata
    with open(os.path.join(base_dir_step1, "stitched_cut_efficient.txt"), "w") as f:
        for file in os.listdir(os.path.join(base_dir_step1, "stitched_cut_efficient")):
            f.write(f"{file}\n")
    with open(os.path.join(base_dir_step1, "final_stitched_cut_efficient.txt"), "w") as f:
        for file in os.listdir(os.path.join(base_dir_step1, "final_stitched_cut_efficient")):
            f.write(f"{file}\n")

    #Step 2: Artefact reject
    print("Artefact rejecting...")
    base_dir_step2 = os.path.join(base_dir,"2_artefact_rejected")
    iterate_artefact_rejection(base_dir_step2)

    #Step 3: Preprocess
    print("Preprocessing...")
    base_dir_step3 = os.path.join(base_dir, "3_preprocessing")
    iterate_preprocess(base_dir_step3) #20 here is a minimum segment size of 20 seconds
    #Metadata
    with open(os.path.join(base_dir_step3, "preprocessed_files.txt"), "w") as f:
        for file in os.listdir(os.path.join(base_dir_step3, r"preprocessed_files")):
            f.write(f"{file}\n")
    
    #Step 4: Variable extraction
    print("Extracting vars...")
    base_dir_step4 = os.path.join(base_dir, r"99_extract_all_variables")
    extract_vars_iterate(base_dir_step4)

    #Step 5: Build GUI to allow user to run stats
    print("Building GUI...")
    stat_run(base_dir)


if __name__ == "__main__":
    setting = "postop" #Can be "intraop" or "postop"
    run_all(setting)

    
