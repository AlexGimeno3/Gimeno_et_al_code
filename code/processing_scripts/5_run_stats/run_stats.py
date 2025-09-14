"""
This is a script that should run a huge set of stats, using the stat_object, which stores all the data.
"""
from stat_gui import stat_gui
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
utils_dir = os.path.join(current_dir, '..', '..', 'utils')
sys.path.append(utils_dir)
from utils import config
vars_dict = config.get_vars_dict()

import time
from datetime import timedelta
import glob

def stat_run(base_dir):
    """
    Here, time_period can be "intraop" or "postop"
    """
    start_time = time.time()
    print("Starting analysis...")
    if vars_dict["windows"] == [["surgery"]]:
        timeseries_file_dir = base_dir + r"\99_extract_all_variables\time_series_data_eeg"
    else:
        timeseries_file_dir = base_dir + r"\99_extract_all_variables\time_series_data_eeg"
    
    if vars_dict["clear"]:
        base_1 = os.path.split(timeseries_file_dir)[0]
        base_cache = os.path.join(base_1, "pkl_cache")
        if not os.path.exists(base_cache):
            os.mkdir(base_cache)
        files = glob.glob(os.path.join(base_cache, "*"))
        for f in files:
            try:
                os.remove(f)
            except Exception as e:
                print(f"Error removing {f}: {e}")

    my_gui = stat_gui(timeseries_file_dir=timeseries_file_dir)
    my_gui.stats_all()

    # Calculate elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time
    # Format the elapsed time
    time_formatted = str(timedelta(seconds=int(elapsed_time)))
    # Print completion message with time
    print(f"Analysis complete! Total time: {time_formatted}")


if __name__ == "__main__":
      #stat_run()
      pass
