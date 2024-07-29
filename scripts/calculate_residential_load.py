import pandas as pd
import numpy as np


if __name__ == "__main__":
    
    # load the number of residential structures
    res_structures = pd.read_csv(snakemake.input.res_structures)
    res_elec_load = pd.read_csv(snakemake.input.res_elec_load)
    res_heat_load = pd.read_csv(snakemake.input.res_heat_load)
    
    total_elec_load = res_elec_load*res_structures.T.loc['n_units']
    total_heat_load = res_heat_load*res_structures.T.loc['n_units']
    
    