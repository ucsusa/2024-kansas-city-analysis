import pandas as pd
import numpy as np


if __name__ == "__main__":
    
    # load price data
    electricity_price = float(snakemake.config['retail_price_elec'])
    gas_price = float(snakemake.config['retail_price_gas'])  # price in $/Mcf
    mcf_to_kwh = 1/297.59333 
    gas_price = gas_price*mcf_to_kwh  # price in kWh
    
    # load the number of residential structures
    res_structures = pd.read_csv(snakemake.input.res_structures, index_col=0)
    res_elec_load = pd.read_csv(snakemake.input.elec_load, parse_dates=True, index_col='timestamp')
    res_heat_load = pd.read_csv(snakemake.input.heat_load, parse_dates=True, index_col='timestamp')
    
    # load energy expenses from the LEAD data
    expenses = pd.read_csv(snakemake.input.res_energy_expenses, index_col='BLD')
    
    
    rescaled_elec_load = (res_elec_load.div(res_elec_load.sum(),axis=1)*
                         (res_elec_load.columns.map(expenses['ELEP*UNITS'])/electricity_price))
    
    rescaled_heat_load = res_heat_load.copy()
    # rescaled_heat_load = (res_heat_load.div(res_heat_load.sum(),axis=1)*
    #                      (res_heat_load.columns.map(expenses['GASP*UNITS'])/gas_price))
    
    total_elec_load = rescaled_elec_load*res_structures.T.loc['n_units']
    total_heat_load = rescaled_heat_load*res_structures.T.loc['n_units']
    
    total_elec_load.to_csv(snakemake.output.rescaled_elec_load)
    total_heat_load.to_csv(snakemake.output.rescaled_heat_load)
    
    
    
    