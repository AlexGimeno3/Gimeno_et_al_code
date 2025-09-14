import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
utils_dir = os.path.join(current_dir, '..', '..', 'utils')
sys.path.append(utils_dir)
spec_path = os.path.join(utils_dir,"NEURAL_py_EEG")
sys.path.insert(0, spec_path)
from utils.processing_parameters_intraop import vars_dict
from spectral_features import get_band_power, spectral_diff, spectral_edge_frequency, spectral_entropy
import neurokit2 as nk
import numpy as np
utils_path = r'D:\utils'
sys.path.append(utils_path)
from EntropyHub import MSEn, MSobject
import numpy as np
from scipy.stats import linregress

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


#1: Spectral power
@add_variable_metadata(variable_name="absolute_spectral_power", units="μV²/Hz", window_size=60, slide_size=30, subfunctions={"delta":default_stats, "theta":default_stats, "alpha":default_stats, "beta":default_stats, "a_d_ratio":default_stats, "t_d_ratio":default_stats})
def spectral_power_abs(data,raw_obj):
    """
    An example function that will be applied to a time series.

    Input:
    - data (arr of float): the array of voltages upon which we would like to apply the function of interest
    - raw_obj (int): alex_raw_efficient object

    Output:
    - delta_power_abs_out (float): absolute delta (0.5-3 Hz) power in our data (units of uV^2/Hz if abs, unitless if relative)
    """
    #window_size of 20 secs chosen after discussion with Dr. Andrezjak, but I will need to double-check. slide_size of 10 chosen since this is an overlap of 50% (which is the default in the NEURAL_py toolbox)
    fs = raw_obj.sampling_rate
    delta = get_band_power(data, fs, freq_bands[0], spectral_power_abs.window_size, relative = False)
    theta = get_band_power(data, fs, freq_bands[1], spectral_power_abs.window_size, relative = False)
    alpha = get_band_power(data, fs, freq_bands[2], spectral_power_abs.window_size, relative = False)
    beta = get_band_power(data, fs, freq_bands[3], spectral_power_abs.window_size, relative = False)
    a_d_ratio = alpha/delta
    t_d_ratio = theta/delta
    return delta, theta, alpha, beta, a_d_ratio, t_d_ratio
    
#3: rEEG
@add_variable_metadata(variable_name="raw_rEEG", units="μV", window_size=2, slide_size=2, subfunctions=["mean", "median", "STDEV", "coefficient_of_variance", "5th_percentile", "95th_percentile", "width", "mann_kendall_tau", "theil_sen_slope", "theil_sen_y_int", "lin_reg", "lin_reg_y"])
def raw_rEEG(data, raw_obj):
    data_point = np.percentile(data, 95) - np.percentile(data, 5)
    return data_point


@add_variable_metadata(variable_name = "raw_rEEG_prop", units = "", window_size = 60, slide_size = 30, subfunctions = {"_0_10": default_stats, "_10_25": default_stats, "_25_50": default_stats, "_50_100": default_stats, "_over_100": default_stats})
def raw_rEEG_props(data,raw_obj):
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
    
#2: Spectral difference
@add_variable_metadata(variable_name="spectral_difference", units="μV²/Hz", window_size=60, slide_size=30, subfunctions={"delta":default_stats, "theta":default_stats, "alpha":default_stats, "beta":default_stats})
def spectral_difference(data, raw_obj):
    #Copied from NEURAL_parameters.py
    params_st = {}
    params_st["L_window"] = 2 #2 sec windows that will give us differences
    params_st["window_type"] = 'hamm' #Hamming window
    params_st["overlap"] = 50 #50% overlap between windows
    params_st["freq_bands"] = freq_bands
    
    delta, theta, alpha, beta = spectral_diff(data, raw_obj.sampling_rate, None, params_st)
    return delta, theta, alpha, beta

#4: Edge frequency
@add_variable_metadata(variable_name="edge_frequency", units="Hz", window_size=60, slide_size=30, subfunctions=default_stats)
def edge_frequency(data, raw_obj):
    params_st = {}
    params_st["L_window"] = 2 #2 sec windows that will give us differences
    params_st["window_type"] = 'hamm' #Hamming window
    params_st["overlap"] = 50 #50% overlap between windows
    params_st["method"] = 'psd' #Default for spectral analysis in E:\EEG_Pipeline_Continuous\99_extract_all_variables\NEURAL_py_EEG\NEURAL_parameters.py
    params_st["SEF"] = 0.95 #Perecentile we want to use to determing SEF
    params_st["total_freq_bands"] = vars_dict["full_range"]
    sef = spectral_edge_frequency(data, raw_obj.sampling_rate, None, params_st)
    return sef

#5: Multiscale entropy:
@add_variable_metadata(variable_name="multiscale_entropy", units="", window_size=60, slide_size=30, subfunctions={"auc": default_stats, "max_value": default_stats, "slope_fine_values": default_stats, "slope_coarse_values": default_stats})
def multiscale_entropy(data, raw_obj):
    # Compute multiscale entropy (MSE)
    scales = list(range(1, 21))  # Scales τ = 1 to 20
    scales = np.array(scales)
    msen_obj = MSobject('SyDyEn')
    mse, x =  MSEn(data, msen_obj, Scales=20, Methodx='coarse', Plotx=False)

    # 1. Area under the MSE curve
    auc = x

    # 2. Maximum value of the MSE curve
    max_value = np.max(mse)

    # 3. Average slope of fine scales (τ = 1 to 5)
    fine_indices = np.where((scales >= 1) & (scales <= 5))[0]
    slope_fine_values, _, _, _, _ = linregress(scales[fine_indices], mse[fine_indices])

    # 4. Average slope of coarse scales (τ = 6 to 20)
    coarse_indices = np.where((scales >= 6) & (scales <= 20))[0]
    slope_coarse_values, _, _, _, _ = linregress(scales[coarse_indices], mse[coarse_indices])

    return auc, max_value, slope_fine_values, slope_coarse_values

#6: Relative spectral power:
@add_variable_metadata(variable_name="relative_spectral_power", units="", window_size=60, slide_size=30, subfunctions={"delta":default_stats, "theta":default_stats, "alpha":default_stats, "beta":default_stats})
def relative_spectral_power(data, raw_obj):
    """
    Calculate relative spectral power in different frequency bands.

    Input:
    - data (arr of float): the array of voltages upon which we would like to apply the function of interest
    - raw_obj (int): alex_raw_efficient object

    Output:
    - tuple of (delta, theta, alpha, beta) powers:
        - if relative=False: absolute power in μV²/Hz
        - if relative=True: relative power (unitless, normalized by total power)
    """
    fs = raw_obj.sampling_rate
    
    delta = get_band_power(data, fs, freq_bands[0], relative_spectral_power.window_size, relative=False)
    theta = get_band_power(data, fs, freq_bands[1], relative_spectral_power.window_size, relative=False)
    alpha = get_band_power(data, fs, freq_bands[2], relative_spectral_power.window_size, relative=False)
    beta = get_band_power(data, fs, freq_bands[3], relative_spectral_power.window_size, relative=False)
    
    total_power = delta + theta + alpha + beta
    return (delta/total_power, theta/total_power, alpha/total_power, beta/total_power)


@add_variable_metadata(variable_name="shannon_entropy", units="units", window_size=60, slide_size=30, subfunctions={"delta":default_stats, "theta":default_stats, "alpha":default_stats, "beta":default_stats})
def shannon_entropy(data,raw_obj):
    feat_name = 42 #Not used
    params_st = {}
    params_st["L_window"] = 2 #2 sec windows that will give us differences
    params_st["window_type"] = 'hamm' #Hamming window
    params_st["overlap"] = 50 #50% overlap between windows
    params_st["freq_bands"] = freq_bands
    params_st["method"] = "PSD"
    entropy_0, entropy_1, entropy_2, entropy_3 = spectral_entropy(data, raw_obj.sampling_rate, feat_name, params_st)
    return (entropy_0, entropy_1, entropy_2, entropy_3)


@add_variable_metadata(variable_name="higuchi_fractal_dimension_raw", units="units", window_size=60, slide_size=30, subfunctions=default_stats)
def higuchi_fractal_dimension_raw(data, raw_obj):
    hfd_value, _ = nk.fractal_higuchi(data, k_max=10)
    return hfd_value







# #8: A statistic (from Adrzejak, 2011)
# @add_variable_metadata(variable_name="nonlinear_a_statistic", units="", window_size=30*60, slide_size=30*60, subfunctions=default_stats)
# def a_statistic(data, raw_obj):
#     # Downsample
#     downsample_factor = 8
#     srate = raw_obj.sampling_rate // downsample_factor
#     data = data[::downsample_factor]
#     reference_srate = 50

#     # Parameters
#     m = 6
#     tau = int(srate / reference_srate * 10)
#     k = 5
#     h = int(srate / reference_srate * 5)
#     w = int(srate / reference_srate * 30)
#     num_surrogates = 19
#     batch_size = 2500

#     # Normalize data on GPU
#     data_gpu = cp.asarray(data, dtype=cp.float32)
#     data_gpu = (data_gpu - cp.mean(data_gpu)) / cp.std(data_gpu)
#     n = len(data_gpu)
#     eta = m * tau

#     # Delay vector creation on GPU
#     indices = cp.arange(eta, n - h, dtype=cp.int32)
#     delay_vectors = cp.stack([data_gpu[indices - (m - i - 1) * tau] for i in range(m)], axis=1)

#     def batch_prediction_error_gpu(delay_vectors, data_gpu, indices, batch_size=batch_size):
#         n_vectors = len(delay_vectors)
#         errors = cp.zeros(n_vectors, dtype=cp.float32)
        
#         for batch_start in tqdm(range(0, n_vectors, batch_size), desc="Processing batches", leave=False):
#             batch_end = min(batch_start + batch_size, n_vectors)
#             batch_vectors = delay_vectors[batch_start:batch_end]
            
#             # Compute pairwise distances
#             distances = cp.linalg.norm(batch_vectors[:, None, :] - delay_vectors[None, :, :], axis=2)
            
#             # Apply proper Theiler window correction
#             for i in tqdm(range(batch_end - batch_start), desc="Applying Theiler Correction", leave=False):
#                 idx = i + batch_start
#                 mask = cp.abs(indices - indices[idx]) < w
#                 distances[i][mask] = cp.inf
                
#             #print("Getting 5 nearest neighbors for each delay vector...")
#             neighbors = cp.argpartition(distances, k, axis=1)[:, :k]
            
#             # Calculate prediction error according to equation 9.12
#             true_horizon = data_gpu[indices[batch_start:batch_end] + h]
            
#             # Progress bar for inner loop
#             for i in tqdm(range(batch_end - batch_start), desc="Computing prediction errors", leave=False):
#                 neighbor_horizons = data_gpu[indices[neighbors[i]] + h]
#                 errors[batch_start + i] = true_horizon[i] - cp.mean(neighbor_horizons)
                
#         return cp.mean(errors ** 2).get()  # Root mean square as per equation 9.13


#     def phase_randomize_gpu():
#         fft_vals = cp.fft.fft(data_gpu)
#         magnitudes = cp.abs(fft_vals)
#         phases = cp.random.uniform(0, 2 * cp.pi, len(data_gpu))
#         random_fft = magnitudes * cp.exp(1j * phases)
#         return cp.real(cp.fft.ifft(random_fft))

#     # Calculate E statistic for original data
#     e_true = batch_prediction_error_gpu(delay_vectors, data_gpu, indices)

#     # Surrogate calculations
#     e_surrogates = []
#     for _ in tqdm(range(num_surrogates)):
#         surrogate_data = phase_randomize_gpu()
#         surrogate_vectors = cp.stack([surrogate_data[indices - (m - i - 1) * tau] for i in range(m)], axis=1)
#         e_surrogates.append(batch_prediction_error_gpu(surrogate_vectors, surrogate_data, indices))

#     return np.mean(e_surrogates) - e_true