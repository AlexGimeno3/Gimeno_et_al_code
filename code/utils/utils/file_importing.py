"""
The point of this module is to contain all the functions I have built or used to import file names
"""
import os
import pickle
import functools
import mne
        
@functools.lru_cache(maxsize=8)
def import_pkl(file_directory=None, file_name=None, full_file_path = None, parts_or_path = "parts", file_path=None):
    r"""
    Loads a previously saved pkl file
    Inputs:
    - file_directory (str): the directory the file can be found in
    - file_name (str): the name of the file. I chose to split these two inputs up in anticipation of looping uses of this function
    - full_file_path (str): the full path of a file (e.g.,FIS104b\FIS104_stitched_filtered.pkl). Accessed when parts_or_path equals path.
    - parts_or_path (str): takes "parts" or "path". If parts, the function will expect that file_directory and file_name were provided; if "path", the function will expect that a file path was given

    Outputs:
    - alex_raw_object (alex_raw object): the alex_raw object contained in the .pkl file
    """
    if parts_or_path == "parts":
        file_path = os.path.join(file_directory, file_name)
    elif parts_or_path == "path":
        file_path = full_file_path

    if not os.path.exists(file_path):
        raise ValueError(f"File {file_path} does not exist.")
    with open(file_path, 'rb') as f:
        print(f"Loading {file_path}")        
        data = pickle.load(f)
        change_Text = f".pkl file loaded from {file_path}"
        data.add_change_log(change_Text)

    if hasattr(data, 'percentiles_envelope'):
        data.percentiles_envelope.parent = data  # Restore the reference to the parent
    data.file_path_current = file_path
    return data

def get_band_obj_from_original(original_obj, band):
    """
    Used in 4_variable_funtions.py to load a signal filtered to a certain frequency band.
    Inputs:
    - original_obj (alex_raw_efficient object)
    - band (str): can be "delta", "theta", "alpha", "beta", or "original" (the unfiltered signal, but in our temporal bound of interest)

    Outputs:
    - band_obj (alex_raw_efficient object): the filtered alex_raw_efficient object of interest
    """
    if band == "":
        return original_obj
    else:
        band_char = band[0]
        base_path = original_obj.file_path_current

        file_base_name = base_path.split('\\')[-1]
        last_pp_idx = file_base_name.rindex('pp')
        dot_idx = file_base_name.rindex('.')
        file_base_name = file_base_name[:last_pp_idx+2] + band_char + file_base_name[last_pp_idx+2:dot_idx] + file_base_name[dot_idx:]

        path_parts = base_path.split('\\')
        new_base = '\\'.join(path_parts[:-2])
        new_path = f"{new_base}\\preprocessed_files_{band}\\{file_base_name}"

        band_obj = import_pkl(parts_or_path="path", full_file_path = new_path)
        return band_obj

def import_edf_from_path(edf_file_path, downsample = None):
    """
    This function returns an MNE raw object from a given EDF file path.

    Inputs:
    - edf_file_path (str): the entire file path for the EDF file of interest
    
    Returns:
    - raw (raw MNE object): the raw EDF MNE file
    """

    raw = mne.io.read_raw_edf(edf_file_path, preload=True, encoding='latin1') #Preload and encoding added to help deal with corrupted files
    raw.original_sampling_rate = int(raw.info["sfreq"])
    if not downsample is None and not(int(raw.info["sfreq"])==downsample):
        print(f"Downsampling to {downsample} from intial sampling rate of {raw.info["sfreq"]}")
        raw.resample(downsample)
    return raw
