_vars_dict = None

def initialize_config(time_frame):
    global _vars_dict
    if time_frame == "intraop":
        from utils.processing_parameters_intraop import vars_dict
    elif time_frame == "postop":
        from utils.processing_parameters_postop import vars_dict
    else:
        raise ValueError("Please use either 'intraop' or 'postop' for the time frame.")
    
    _vars_dict = vars_dict

def get_vars_dict():
    if _vars_dict is None:
        raise RuntimeError("Config not initialized. Call initialize_config() first.")
    return _vars_dict

# Convenience accessors
def get_windows():
    return get_vars_dict()["windows"]

def get_base_dir():
    return get_vars_dict()["base_dir"]