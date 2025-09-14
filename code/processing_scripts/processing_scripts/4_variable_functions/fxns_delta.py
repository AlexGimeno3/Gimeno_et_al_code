import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
utils_dir = os.path.join(current_dir, '..', '..', 'utils')
sys.path.append(utils_dir)
spec_path = os.path.join(utils_dir,"NEURAL_py_EEG")
sys.path.insert(0, spec_path)
from utils.processing_parameters_intraop import vars_dict

import numpy as np
import neurokit2 as nk
from scipy import signal
from scipy.stats import kurtosis, skew

freq_bands = [vars_dict["delta_range"], vars_dict["theta_range"], vars_dict["alpha_range"], vars_dict["beta_range"]]
default_stats = ["median", "STDEV", "IQR", "mann_kendall_tau", "theil_sen_slope", "theil_sen_y_int", "lin_reg", "lin_reg_y", "hurst"]

#The subfunctions variable in add_variable_metadata defaults to none, but can be edited when defining the metadata to contain every subfunction I want to analyze for that given complete function time series
def add_variable_metadata(variable_name, units, window_size, slide_size, subfunctions=None):
    """Decorator to add metadata about the variable being analyzed. Does NOT get read in extract_variables_segmented"""
    def decorator(func):
        # Add metadata as attributes to the function
        func.variable_name = variable_name
        func.units = units
        func.window_size = window_size #In seconds
        func.slide_size = slide_size #In seconds
        func.subfunctions = subfunctions or []  # Default to an empty list if not provided
        return func
    return decorator

#rEEG
@add_variable_metadata(variable_name="delta_band_rEEG", units="μV", window_size=2, slide_size=2, subfunctions=["mean", "median", "STDEV", "coefficient_of_variance", "5th_percentile", "95th_percentile", "width", "mann_kendall_tau", "theil_sen_slope", "theil_sen_y_int", "lin_reg", "lin_reg_y"])
def rEEG_delta_band(data,raw_obj):
    data_point = np.percentile(data, 95) - np.percentile(data, 5)
    return data_point

@add_variable_metadata(variable_name = "delta_rEEG_prop", units = "", window_size = 60, slide_size = 30, subfunctions = {"_0_10": default_stats, "_10_25": default_stats, "_25_50": default_stats, "_50_100": default_stats, "_over_100": default_stats})
def delta_rEEG_props(data,raw_obj):
    rEEG_window = []
    num_windows = int(60/2) #U60 from window_size, 2 because we take samples every 2 seconds
    mini_window_size = 2*raw_obj.sampling_rate
    #Calculate our rEEG segment
    for i in range(num_windows):
            start_idx = i * mini_window_size
            end_idx = start_idx + mini_window_size
            data_mini = data[start_idx:end_idx]
            window_range = np.percentile(data_mini, 95) - np.percentile(data_mini, 5)
            rEEG_window.append(window_range)
    rEEG_window = np.array(rEEG_window)

    _0_10 = np.sum((rEEG_window >= 0) & (rEEG_window < 10)) / len(rEEG_window)
    _10_25 = np.sum((rEEG_window >= 10) & (rEEG_window < 25)) / len(rEEG_window)
    _25_50 = np.sum((rEEG_window >= 25) & (rEEG_window < 50)) / len(rEEG_window)
    _50_100 = np.sum((rEEG_window >= 50) & (rEEG_window <= 100)) / len(rEEG_window)
    _over_100 = np.sum((rEEG_window > 100)) / len(rEEG_window)

    return _0_10, _10_25, _25_50, _50_100, _over_100


#3: Envelope
@add_variable_metadata(variable_name="delta_envelope", units="μV²", window_size=60, slide_size=30, subfunctions={"delta_env_mean": default_stats, "delta_env_STDEV": default_stats})
def delta_envelope(data, raw_obj):
    env = np.abs(signal.hilbert(data)) ** 2
    delta_env_mean = np.nanmean(env)
    delta_env_STDEV = np.nanstd(env, ddof=1)
    return delta_env_mean, delta_env_STDEV

#5: Power
@add_variable_metadata(variable_name="delta_power", units="μV²", window_size=60, slide_size=30, subfunctions={"mean_power_delta": default_stats, "stdev_delta": default_stats})
def delta_power(data, raw_obj):
    mean_power_delta =  np.nanmean(np.abs(data) ** 2)
    stdev_delta = np.std(data, ddof=1)
    return mean_power_delta, stdev_delta

#6: Skew
@add_variable_metadata(variable_name="delta_skew", units="units", window_size=60, slide_size=30, subfunctions=default_stats)
def delta_skew(data, raw_obj):
    return np.abs(skew(data, nan_policy="omit"))

#4: Kurtosis
@add_variable_metadata(variable_name="delta_kurtosis", units="units", window_size=60, slide_size=30, subfunctions=default_stats)
def delta_kurtosis(data, raw_obj):
    return kurtosis(data, fisher=False, nan_policy='omit')

@add_variable_metadata(variable_name="higuchi_fractal_dimension_delta", units="units", window_size=60, slide_size=30, subfunctions=default_stats)
def higuchi_fractal_dimension_delta(data, raw_obj):
    hfd_value, _ = nk.fractal_higuchi(data, k_max=10)
    return hfd_value

