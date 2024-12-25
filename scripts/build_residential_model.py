import numpy as np
import pandas as pd
import pypsa
import matplotlib.pyplot as plt
from unyt import MWh, kWh

def annuity(r, n):
    return r / (1 - 1 / (1 + r)**n)


if __name__ == "__main__":
    ###############################################################
    # create the network container and snapshots
    ###############################################################
    n = pypsa.Network(name=snakemake.config['community_name'])

    N_days=365
    N_hours=24

    index = pd.date_range(start="2018-01-01", 
                        periods=N_days*N_hours, 
                        freq='h')

    n.set_snapshots(index)

    ###############################################################
    # add energy carriers
    ###############################################################
    n.add(class_name="Carrier", name="grid")
    n.add(class_name="Carrier", name="solar")
    n.add(class_name="Carrier", name="battery")
    n.add(class_name="Carrier", name='net metering')




    ###############################################################
    # add buses
    ###############################################################
    bus_names = ['Residential','CommunitySolar']
    for bus in bus_names:
        n.add(class_name="Bus",
                name=bus,
                carrier='AC')
        

    ###############################################################
    # add link connecting community solar to residential
    ###############################################################
    delivery_price = float(snakemake.config['delivery_price_elec'])

    n.add(class_name='Link',
      name=f"CommSol-Res",
      bus0="CommunitySolar",
      bus1="Residential",
      carrier="AC",
      p_nom_extendable=True,
      marginal_cost=delivery_price,
      )


    ###############################################################
    # add electricity demand
    ###############################################################
    load = pd.read_csv(snakemake.input.rescaled_elec_load, 
                       parse_dates=True, 
                       index_col='timestamp')
    
    load_resampled = load.loc['2018'].resample('h').mean().sum(axis=1)
    load_resampled = load_resampled / 1e3  # kW --> MW


    bus_name = 'Residential'
    n.add(
        class_name="Load",
        name=f"Load {bus_name}",
        bus=bus_name,
        p_set=load_resampled
        )
    
    ###############################################################
    # load weather data
    ###############################################################

    weather = pd.read_csv(snakemake.input.weather, 
                          parse_dates=True, 
                          index_col=0)
    
    ghi = weather['ghi'] / weather['ghi'].max()

    ###############################################################
    # load cost data
    ###############################################################

    costs = pd.read_csv(snakemake.input.costs, 
                        index_col='technology')
    costs *= 1e3  # convert /kW to /MW

    utility_costs = pd.read_csv(snakemake.input.utility_costs, 
                                index_col='technology')
    utility_costs *= 1e3 # convert /kW to /MW

    annuity_adj = annuity(0.07, 20)

    costs = costs.assign(annualized_cost = (costs['OCC']*annuity_adj 
                                            + costs['Fixed O&M']))
    utility_costs = utility_costs.assign(
        annualized_cost = (utility_costs['OCC']*annuity_adj 
                            + utility_costs['Fixed O&M'])
        )
    
    retail_price = snakemake.config['retail_price_elec']
    
    ###############################################################
    # add generators
    ###############################################################

    # residential generators
    for generator in costs.index:
        if generator == 'DistributedWind':
            pass
        else:
            print(generator)
            annualized_cost = costs.at[generator, 'annualized_cost']
            print(annualized_cost)
            
            if generator=='ResPV':
                n.add(class_name='Generator',
                        name=generator,
                        bus=bus_name,
                        carrier="solar",
                        capital_cost=annualized_cost,  # $/kW
                        p_min_pu=ghi,
                        p_max_pu=ghi,
                        p_nom_extendable=True,
                        p_nom_max = 3.0,
                        )
            elif generator=='Residential Battery Storage':
                pass
                n.add(class_name="StorageUnit",
                        name=generator,
                        bus=bus_name,
                        carrier="battery",
                        capital_cost=annualized_cost,  # $/kW
                        p_nom_extendable=True,
                        max_hours=2.5,
                        cyclic_state_of_charge=False,
                        )
        
    # community solar generator
    n.add(class_name="Generator",
        name="Community Solar",
        bus="CommunitySolar",
        carrier="solar",
        capital_cost=utility_costs.at['UtilityPV','annualized_cost'],
        p_min_pu=ghi,
        p_max_pu=ghi,
        p_nom_extendable=True,
        p_nom_max=3.0
        )
    
    # net metering
    for bus in bus_names:
      n.add(class_name="Generator",
            name=f"Net metering {bus}",
            bus=bus,
            carrier='net metering',
            p_min_pu=-1,
            p_max_pu=0.0,
            marginal_cost=retail_price*0.0,
            capital_cost=0.0,
            p_nom_extendable=True)
    
    # electricity imports
    n.add(class_name='Generator',
      name='Evergy Import',
      bus=bus_name,
      carrier='grid',
      capital_cost=0,
      marginal_cost=retail_price,
      p_nom_extendable=True,
      )
        

    n.export_to_netcdf(snakemake.output.residential_model)