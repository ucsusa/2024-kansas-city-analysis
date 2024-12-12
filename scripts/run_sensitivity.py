import pandas as pd
import numpy as np
import pypsa
import matplotlib.pyplot as plt


if __name__ == "__main__":
    solver = snakemake.config['solver']

    n = pypsa.Network(snakemake.input.residential_model)

    try:
        n.optimize(solver_name=solver)
    except:
        print(f'No {solver} solver found. Reverting to HiGHS')
        n.optimize(solver_name='highs')


    data = {# sensitivities/cases
            'rooftop_discount':[],
            'community_solar_discount':[],
            'battery_discount':[],
            'percent_retail_price':[],
            'community_solar_potential':[],
            # outputs
            'rooftop_capacity':[],
            'community_solar_capacity':[],
            'residential_battery_capacity':[],
            'objective_value':[],
            }
    

    discounts = np.arange(0, 1+delta, delta)
    retail_prices = np.linspace(0, 1, 5)
    for discount in discounts:
        for pct_retail in retail_prices:
            n.generators.loc['ResPV', 'capital_cost'] = costs.at['ResPV','annualized_cost'] * (1-discount)
            n.generators.loc['Net metering Residential', 'marginal_cost'] = retail_price*pct_retail
            try:
                n.optimize(solver_name=solver)
            except:
                print(f'No {solver} solver found. Reverting to HiGHS')
                n.optimize(solver_name='highs')
            
            data['discount'].append(discount)
            data['solar_capacity'].append(np.abs(n.generators.p_nom_opt['ResPV']))
            data['battery_capacity'].append(np.abs(n.storage_units.p_nom_opt['Residential Battery Storage']))
            data['percent_retail_price'].append(pct_retail)
            data['objective_value'].append(n.objective)