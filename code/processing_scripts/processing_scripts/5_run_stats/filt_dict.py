from copy import deepcopy
import numpy as np

class filt_dict:
    def __init__(self, fis_obj_dict):
        dict_size = len(list(fis_obj_dict.values()))
        filters_dict = {
            "fis":np.empty(dict_size),
            "ND_filter":np.empty(dict_size),
            "ecc_filter":np.empty(dict_size),
            "machine_type_filter":np.empty(dict_size),
            "ignore_filter":np.empty(dict_size),
            "seizure_filter":np.empty(dict_size),
            "artefact_filter":np.empty(dict_size)
        }
        
        i = 0
        for fis_obj in fis_obj_dict.values():
            filters_dict["fis"][i] = deepcopy(fis_obj.fis)
            filters_dict["ND_filter"][i] = deepcopy(fis_obj.ND_filter)
            filters_dict["ecc_filter"][i] = deepcopy(fis_obj.ecc_filter)
            filters_dict["machine_type_filter"][i] = deepcopy(fis_obj.machine_type_filter)
            filters_dict["ignore_filter"][i] = deepcopy(fis_obj.ignore_filter)
            filters_dict["seizure_filter"][i] = deepcopy(fis_obj.seizure_filter)
            filters_dict["artefact_filter"][i] = deepcopy(fis_obj.artefact_filter)
            i=i+1
        
        self.filters_dict = filters_dict
    
    def get_filter(self, fis):
        """
        Returns the current filter value of a given fis (int).
        """
        i = np.where(self.filters_dict["fis"]==fis)
        return self.filters_dict["ND_filter"][i]==1 and self.filters_dict["ecc_filter"][i]==1 and self.filters_dict["machine_type_filter"][i]==1 and self.filters_dict["ignore_filter"][i]==1 and self.filters_dict["seizure_filter"][i]==1 and self.filters_dict["artefact_filter"][i]==1
    
    def set_filter(self, fis, filter_name, setting):
        """
        Allows us to set a filter value.
        Inputs:
        - fis (int): fis whose filter we would like to set
        - filter_name (str): string of the filter we want to set. Options are "ND", "ecc", "machine_type", "ignore", "seizure", "artefact"
        - setting (int): 0 or 1 (0 to not include in analysis, 1 to include in analysis)
        """
        filter_string = filter_name+"_filter"
        i = np.where(self.filters_dict["fis"]==fis)
        self.filters_dict[filter_string][i]=setting

