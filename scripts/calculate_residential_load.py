import pandas as pd
import numpy as np


if __name__ == "__main__":
    
    # load the number of residential structures
    res_structures = pd.read_csv(snakemake.input.res_structures, index_col=0)
    res_elec_load = pd.read_csv(snakemake.input.res_elec_load, parse_dates=True, index_col='timestamp')
    res_heat_load = pd.read_csv(snakemake.input.res_heat_load, parse_dates=True, index_col='timestamp')
    
    total_elec_load = res_elec_load*res_structures.T.loc['n_units']
    total_heat_load = res_heat_load*res_structures.T.loc['n_units']
    
    