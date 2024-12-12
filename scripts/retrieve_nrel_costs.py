import numpy as np
import pandas as pd
from nrelpy.atb import ATBe

if __name__ == "__main__":
    atb_params = snakemake.config['atb_params']

    atb = ATBe(atb_params['atb_year'])
    df = atb.raw_dataframe
    new_selection = df[
                        (df['core_metric_case']==atb_params['case'])
                        &(df['scale']==atb_params['scale'])
                        &(df['maturity']==atb_params['maturity'])
                        &(df['scenario']==atb_params['scenario'])
                        &(df['core_metric_variable']==atb_params['cost_year'])
                        &(df['default']==1)
                        &(df['crpyears']==atb_params['crp'])
                        &(df['core_metric_parameter'].isin(['Fixed O&M', 'OCC']))
                        ]

    pivot_table = new_selection.pivot_table(index=['technology'],
                                            columns=['core_metric_parameter'],
                                            values='value')
    
    pivot_table.to_csv(snakemake.output.costs)