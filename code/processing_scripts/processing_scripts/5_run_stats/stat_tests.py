from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
from numpy.dtypes import StringDType


def var_dict():
    my_general_dict = {
        "absolute_spectral_power_a_d_ratio": "alpha-delta ratio",
        "absolute_spectral_power_t_d_ratio":"theta-delta ratio",
        "absolute_spectral_power_alpha": "absolute alpha power",
        "absolute_spectral_power_beta":"absolute beta power",
        "absolute_spectral_power_delta":"absolute delta power",
        "absolute_spectral_power_theta":"absolute theta power",
        "relative_spectral_power_alpha": "relative alpha power",
        "relative_spectral_power_beta":"relative beta power",
        "relative_spectral_power_delta":"relative delta power",
        "relative_spectral_power_theta":"relative theta power",
        "edge_frequency":"edge frequency (at 95%)",
        "raw_rEEG":"full-spectrum rEEG signal",
        "higuchi_fractal_dimension_alpha":"alpha band higuchi fractal dimension",
        "higuchi_fractal_dimension_beta":"beta band higuchi fractal dimension",
        "higuchi_fractal_dimension_theta":"theta band higuchi fractal dimension",
        "higuchi_fractal_dimension_delta":"delta band higuchi fractal dimension",
        "higuchi_fractal_dimension_raw":"full-spectrum higuchi fractal dimension",
        "multiscale_entropy_auc":"area under the curve for multiscale entropy",
        "multiscale_entropy_max_value":"multiscale entropy max value",
        "multiscale_entropy_slope_coarse_values":"multiscale entropy slope for coarse values",
        "multiscale_entropy_slope_fine_values":"multiscale entropy slope for fine values",
        "raw_rEEG_prop__0_10":"full-spectrum rEEG proportion between 0 and 10 uV",
        "raw_rEEG_prop__10_25":"full-spectrum rEEG proportion between 10 and 25 uV",
        "raw_rEEG_prop__25_50":"full-spectrum rEEG proportion between 25 and 50 uV",
        "raw_rEEG_prop__50_100":"full-spectrum rEEG proportion between 50 and 100 uV",
        "raw_rEEG_prop__over_100":"full-spectrum rEEG proportion over 100 uV",
        "shannon_entropy_alpha":"alpha band Shannon entropy",
        "shannon_entropy_beta":"beta band Shannon entropy",
        "shannon_entropy_delta":"delta band Shannon entropy",
        "shannon_entropy_theta":"theta band Shannon entropy",
        "spectral_difference_alpha":"alpha band spectral difference",
        "spectral_difference_beta":"beta band spectral difference",
        "spectral_difference_delta":"delta band spectral difference",
        "spectral_difference_theta":"theta band spectral difference",
        }
    my_alpha_dict = {
        "alpha_band_rEEG":"alpha band rEEG",
        "alpha_envelope_alpha_env_mean": "alpha band envelope mean value",
        "alpha_envelope_alpha_env_STDEV":"alpha band envelope standard deviation",
        "alpha_kurtosis":"alpha band kurtosis",
        "alpha_power_mean_power_alpha":"mean alpha band power",
        "alpha_power_stdev_alpha":"standard deviation of alpha band power",
        "alpha_rEEG_prop__0_10":"alpha band rEEG proportion between 0 and 10 uV",
        "alpha_rEEG_prop__10_25":"alpha band rEEG proportion between 10 and 25 uV",
        "alpha_rEEG_prop__25_50":"alpha band rEEG proportion between 25 and 50 uV",
        "alpha_rEEG_prop__50_100":"alpha band rEEG proportion between 50 and 100 uV",
        "alpha_rEEG_prop__over_100":"alpha band rEEG proportion over 100 uV",
        "alpha_skew":"alpha skew"
    }
    my_beta_dict = {
        "beta_band_rEEG":"beta band rEEG",
        "beta_envelope_beta_env_mean": "beta band envelope mean value",
        "beta_envelope_beta_env_STDEV":"beta band envelope standard deviation",
        "beta_kurtosis":"beta band kurtosis",
        "beta_power_mean_power_beta":"mean beta band power",
        "beta_power_stdev_beta":"standard deviation of beta band power",
        "beta_rEEG_prop__0_10":"beta band rEEG proportion between 0 and 10 uV",
        "beta_rEEG_prop__10_25":"beta band rEEG proportion between 10 and 25 uV",
        "beta_rEEG_prop__25_50":"beta band rEEG proportion between 25 and 50 uV",
        "beta_rEEG_prop__50_100":"beta band rEEG proportion between 50 and 100 uV",
        "beta_rEEG_prop__over_100":"beta band rEEG proportion over 100 uV",
        "beta_skew":"beta skew"
    }
    my_theta_dict = {
        "theta_band_rEEG":"theta band rEEG",
        "theta_envelope_theta_env_mean": "theta band envelope mean value",
        "theta_envelope_theta_env_STDEV":"theta band envelope standard deviation",
        "theta_kurtosis":"theta band kurtosis",
        "theta_power_mean_power_theta":"mean theta band power",
        "theta_power_stdev_theta":"standard deviation of theta band power",
        "theta_rEEG_prop__0_10":"theta band rEEG proportion between 0 and 10 uV",
        "theta_rEEG_prop__10_25":"theta band rEEG proportion between 10 and 25 uV",
        "theta_rEEG_prop__25_50":"theta band rEEG proportion between 25 and 50 uV",
        "theta_rEEG_prop__50_100":"theta band rEEG proportion between 50 and 100 uV",
        "theta_rEEG_prop__over_100":"theta band rEEG proportion over 100 uV",
        "theta_skew":"theta skew"
    }
    my_delta_dict = {
        "delta_band_rEEG":"delta band rEEG",
        "delta_envelope_delta_env_mean": "delta band envelope mean value",
        "delta_envelope_delta_env_STDEV":"delta band envelope standard deviation",
        "delta_kurtosis":"delta band kurtosis",
        "delta_power_mean_power_delta":"mean delta band power",
        "delta_power_stdev_delta":"standard deviation of delta band power",
        "delta_rEEG_prop__0_10":"delta band rEEG proportion between 0 and 10 uV",
        "delta_rEEG_prop__10_25":"delta band rEEG proportion between 10 and 25 uV",
        "delta_rEEG_prop__25_50":"delta band rEEG proportion between 25 and 50 uV",
        "delta_rEEG_prop__50_100":"delta band rEEG proportion between 50 and 100 uV",
        "delta_rEEG_prop__over_100":"delta band rEEG proportion over 100 uV",
        "delta_skew":"delta skew"
    }
    merged_dict = {**my_general_dict, **my_delta_dict,  **my_theta_dict,  **my_alpha_dict, **my_beta_dict}
    return merged_dict

def stat_dict():
    my_dict = {
        "mean":"mean",
        "median":"median",
        "sem":"standard error of the mean (SEM)",
        "stdev" : "STDEV",
        "iqr" : "IQR",
        "perc_5":"5th percentile value",
        "perc_95":"95th percentile value",
        "lin_reg_slope":"linear regression slope",
        "lin_reg_yint":"linear regression y-intercept",
        "lin_reg_rmse":"RMSE for linear regression",
        "theil_sen_slope":"Theil-Sen slope",
        "theil_sen_yint":"Theil-Sen y-intercept",
        "theil_sen_rmse":"RMSE for Theil-Sen line of best fit",
        "mann_kendall_tau":"Mann-Kendall tau value",
        "mann_kendall_p":"Mann-Kendall p-value",
        "timeseries_num_points":"Number of data points",
        "timeseries_duration_mins":"Last minute - first minute"
    }
    return my_dict

def param_tests():
    epsilon = np.finfo(float).eps

    tests = {
        "mean": lambda x, y: np.mean(x) if len(x) > 0 else np.nan,
        "median": lambda x, y: np.median(x) if len(x) > 0 else np.nan,
        "sem": lambda x, y: np.std(x, ddof=1) / np.sqrt(len(x)) if len(x) > 1 else 0,
        "stdev": lambda x, y: np.std(x, ddof=1) if len(x) > 1 else 0,
        "iqr": lambda x, y: np.subtract(*np.percentile(x, [75, 25])) if len(x) > 0 else 0,
        "perc_5": lambda x, y: np.percentile(x, 5) if len(x) > 0 else np.nan,
        "perc_95": lambda x, y: np.percentile(x, 95) if len(x) > 0 else np.nan,

        "lin_reg_slope": lambda x, y: np.polyfit(x, y/3600, 1)[0] if len(x) >= 2 and len(np.unique(x)) > 1 else 0,
        "lin_reg_yint": lambda x, y: np.polyfit(x, y/3600, 1)[1] if len(x) >= 2 and len(np.unique(x)) > 1 else np.mean(y) if len(y) > 0 else np.nan,
        "lin_reg_rmse": lambda x, y: np.sqrt(np.mean((y/3600 - np.polyval(np.polyfit(x, y/3600, 1), x))**2)) if len(x) >= 2 and len(np.unique(x)) > 1 else 0, 

        "theil_sen_slope": lambda data, time: stats.theilslopes(data, time/3600)[0] if len(np.unique(data)) > 1 else 0,
        "theil_sen_yint": lambda data, time: stats.theilslopes(data, time/3600)[1] if len(np.unique(data)) > 1 else np.mean(data),
        "theil_sen_rmse": lambda data, time: np.sqrt(np.mean((data - (stats.theilslopes(data, time/3600)[0] * time/3600 + stats.theilslopes(data, time/3600)[1]))**2)) if len(np.unique(data)) > 1 else 0,
        "mann_kendall_tau": lambda data, time: stats.kendalltau(time/3600, data)[0] if len(np.unique(data)) > 1 else 0,
        "mann_kendall_p" : lambda data, time: stats.kendalltau(time/3600, data)[1] if len(np.unique(data)) > 1 else 0,
        "timeseries_num_points" : lambda data, time: len(data),
        "timeseries_duration_mins" : lambda data, time: (time[-1]-time[0])/60
    }
    
    return tests

def create_boxplot(title, nd_data, non_nd_data, p_value):
    # Create a new figure window for each boxplot
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    # Create boxplot
    boxplot_data = [nd_data, non_nd_data]
    box = ax.boxplot(boxplot_data, patch_artist=True, notch=True, labels=['ND', 'Non-ND'])
    # Add individual data points (jitter)
    for i, data in enumerate([nd_data, non_nd_data]):
        # Add jitter for better visualization
        x = np.random.normal(i+1, 0.04, size=len(data))
        ax.plot(x, data, 'r.', beta=0.3, markersize=4)
    # Customize boxplot colors
    colors = ['lightblue', 'lightgreen']
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
    # Add title with p-value
    p_value_str = f"{p_value:.4f}" if p_value is not None else "N/A"
    plt.title(f"{title} by ND Status\np = {p_value_str}", fontsize=12)
    plt.ylabel(title)
    # Add sample size information
    plt.annotate(f"n = {len(nd_data)}", xy=(1, ax.get_ylim()[0]), xycoords='data', ha='center')
    plt.annotate(f"n = {len(non_nd_data)}", xy=(2, ax.get_ylim()[0]), xycoords='data', ha='center')
    # Adjust layout
    plt.tight_layout()
    # Show the figure
    plt.show(block=False)

def change_resolution(arr, digits = 3):
    """ 
    Function that takes in an array of values and returns them to a given number of decimal places with scientific notation.

    Inputs:
    - arr (arr of int or float): array containing the values to be transformed
    - number of digits of significance we want
    
    Outputs:
    - out_arr (arr of str): array containing a string of the modified values
    """
    out_arr = np.empty(len(arr), dtype=StringDType())
    for i in range(len(arr)):
        val = arr[i]
        if val == 0:
            suffix = "0"*digits
            out_arr[i] = f"0.{suffix}"
        elif abs(val)<1/(10**digits):
            out_arr[i] = str(np.format_float_scientific(val, digits+1))
        else:
            out_arr[i] = f"{val:.3f}"
    
    return tuple(out_arr.tolist())

def alex_mwu(group_1_vals, group_2_vals, boxplot = False, var_name = "", param = ""):
    """
    var_name example: delta power
    param example: standard deviation
    """
    # Non-parametric test (Mann-Whitney U)
    u_stat, p_value = stats.mannwhitneyu(group_1_vals, group_2_vals, alternative='two-sided', method="auto")
    g1_med = np.median(group_1_vals)
    g1_iqr = np.percentile(group_1_vals, 0.75)- np.percentile(group_1_vals, 0.25)
    g2_med = np.median(group_2_vals)
    g2_iqr = np.percentile(group_2_vals, 0.75)- np.percentile(group_2_vals, 0.25)
    if boxplot:
        title = f"MWU test run on variable {var_name} {param}"
        create_boxplot(title, group_1_vals, group_2_vals, p_value)
    return g1_med, g1_iqr, g2_med, g2_iqr, u_stat, p_value, None

def alex_students(group_1_vals, group_2_vals, boxplot = False, var_name = "", param = ""):
    # Check variance equality
    _, p_var = stats.levene(group_1_vals, group_2_vals)
    equal_var = p_var >= 0.05
    g1_mean = np.mean(group_1_vals)
    g1_sem = np.std(group_1_vals, ddof=1) / np.sqrt(len(group_1_vals))
    g2_mean = np.mean(group_2_vals)
    g2_sem = np.std(group_2_vals, ddof=1) / np.sqrt(len(group_2_vals))
    
    # Parametric test (t-test)
    t_stat, p_value = stats.ttest_ind(group_1_vals, group_2_vals, equal_var=equal_var)
    if boxplot:
        if equal_var:
            title = f"Student's t-test run on variable {var_name} {param}, equal variances assumed"
        else:
            title = f"Student's t-test run on variable {var_name} {param}, equal variances not assumed"
        create_boxplot(title, group_1_vals, group_2_vals, p_value)
    
    return g1_mean, g1_sem, g2_mean, g2_sem, t_stat, p_value, equal_var