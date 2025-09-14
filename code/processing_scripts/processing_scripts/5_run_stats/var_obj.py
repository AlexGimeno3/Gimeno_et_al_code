import numpy as np

class var_object():
    def __init__(self, name, vals, times, srate):
        """
        A specialized object that stores information about each variable of interest.

        Inputs:
        - name (str): the name of the variable
        - vals (float arr): all the values over time for that given variable
        - times (float arr): all the time values for this variable in SECONDS from surgery end
        - srate (float): the sampling rate of this variable
        """
        
        self.name = name
        self.values = np.array(vals, dtype=np.float32)
        self.times = np.array(times, dtype=np.float32)
        self.srate = srate
    
    def get_vals(self, start_time, end_time, op_status):
        """
        This function returns all values and associated times between start_time and end_time
        Inputs:
        - start_time (float): time from end of surgery where we would like to start our data retrieval IN SECONDS
        - end_time (float): time from end of surgery where we would like to end our data retrieval IN SECONDS
        - op_status: can be "intraop" or "postop"
        """
        # Find indices where times are between start_time and end_time
        if op_status == "postop":
            start_time = start_time*3600
            end_time = end_time*3600
        time_mask = (self.times >= start_time) & (self.times <= end_time)
        
        # Get the times and values using the mask
        ret_times = self.times[time_mask]
        ret_vals = self.values[time_mask]
        
        return ret_times, ret_vals
