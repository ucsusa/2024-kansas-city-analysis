import pandas as pd
import numpy as np
import pypsa
import itertools as it
from tqdm import tqdm
import time

def annuity(r, n):
    return r / (1 - 1 / (1 + r)**n)


if __name__ == "__main__":
    solver = snakemake.config['solver']

    n = pypsa.Network(snakemake.input.residential_model)

    costs = pd.read_csv(snakemake.input.costs, index_col='technology')
    costs *= 1e3  # convert /kW to /MW
    

    utility_costs = pd.read_csv(snakemake.input.utility_costs, index_col='technology')
    utility_costs *= 1e3


    annuity_adj = annuity(float(snakemake.config['interest_rate']), 
                          int(snakemake.config['payback_period']))
    costs = costs.assign(annualized_cost = costs['OCC']*annuity_adj + costs['Fixed O&M'])
    utility_costs = utility_costs.assign(annualized_cost = utility_costs['OCC']*annuity_adj + utility_costs['Fixed O&M'])

    try:
        n.optimize(solver_name=solver)
    except:
        print(f'No {solver} solver found. Reverting to HiGHS')
        n.optimize(solver_name='highs')


    data = {# sensitivities/cases
            'rooftop_discount':[],
            'community_solar_discount':[],
            'percent_retail_price_roof':[],
            'percent_retail_price_comm':[],
            'community_solar_potential':[],
            # outputs
            'rooftop_capacity':[],
            'community_solar_capacity':[],
            'solar_penetration':[],
            'objective_value':[],
            'lcoe':[],
            'savings':[]
            }

    retail_price = snakemake.config['retail_price_elec'] * 1e3

    delta = 0.1
    discounts = np.arange(0, 1+delta, delta)
    retail_prices = np.linspace(0, 1, 9)
    commsol_net_meter = retail_prices.copy()
    community_solar_potential = np.arange(0,3.25,0.25)

    # test ranges
    # discounts = [0, 0.5, 1]
    # retail_prices = np.array([0])
    # commsol_net_meter = [0]
    # community_solar_potential = [0,3]

    sensitivities = list(it.product(discounts, 
                    discounts, 
                    retail_prices, 
                    commsol_net_meter, 
                    community_solar_potential, 
                    ))

    i = 0
    for roof_discount, comm_discount, pct_retail_roof, pct_retail_comm, comm_potential in tqdm(sensitivities, 
                                                                                               leave=True):
        # roof discount
        n.generators.loc['ResPV', 'capital_cost'] = costs.at['ResPV',
                                                             'annualized_cost'] * (1-roof_discount)
        # community solar discount
        n.generators.loc['Community Solar', 'capital_cost'] = utility_costs.at['UtilityPV',
                                                                               'annualized_cost'] * (1-comm_discount)
        # net metering residential
        n.generators.loc['Net metering Residential', 'marginal_cost'] = retail_price*pct_retail_roof
        # net metering commercial
        n.generators.loc['Net metering CommunitySolar', 'marginal_cost'] = retail_price*pct_retail_comm
        # limit community solar capacity
        n.generators.at['Community Solar', 'p_nom_max'] = comm_potential
        try:
            n.optimize(solver_name=solver)
        except:
            print(f'No {solver} solver found. Reverting to HiGHS')
            n.optimize(solver_name='highs')
        
        roof_cap = np.abs(n.generators.p_nom_opt['ResPV'])
        comm_cap = np.abs(n.generators.p_nom_opt['Community Solar'])
        data['rooftop_discount'].append(roof_discount)
        data['community_solar_discount'].append(comm_discount)
        data['percent_retail_price_roof'].append(pct_retail_roof)
        data['percent_retail_price_comm'].append(pct_retail_comm)
        data['community_solar_potential'].append(comm_potential)
        data['rooftop_capacity'].append(roof_cap)
        data['community_solar_capacity'].append(comm_cap)
        data['solar_penetration'].append((roof_cap+comm_cap)/3.0)
        data['objective_value'].append(n.objective)
        lcoe = np.round(n.objective / n.loads_t.p_set.sum().values[0],5)
        data['lcoe'].append(lcoe)
        savings = np.round(1 - lcoe/retail_price,5)
        data['savings'].append(savings)


        timestr = time.strftime("%Y%m%d-%H%M%S")
        if i % 10e3 == 0:
            results_df = pd.DataFrame(data)
            results_df.to_csv(f"results/partial_results/{timestr}_results.csv")

        i += 1
        
        # n.export_to_netcdf(f"data/sensitivity/{timestr}-armourdale-model.nc")

    results_df = pd.DataFrame(data)
    results_df.to_csv(snakemake.output.sensitivity_results)