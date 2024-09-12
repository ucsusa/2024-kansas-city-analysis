import pandas as pd


if __name__ == "__main__":
    
    lead_df = pd.read_csv(snakemake.input.lead_community)
    
    bldg_types = {'1 ATTACHED':"single-family_attached", 
                '1 DETACHED':"single-family_detached",
                '2 UNIT':"multi-family_with_2_-_4_units",
                '3-4 UNIT':"multi-family_with_2_-_4_units",
                '50+ UNIT':"multi-family_with_5plus_units",
                'MOBILE_TRAILER':"mobile_home"}
    
    lead_df = lead_df.replace(bldg_types)
    
    by_unit = lead_df[['BLD','ELEP*UNITS','GASP*UNITS','HINCP*UNITS','UNITS']].groupby(['BLD']).sum()
    
    by_unit = by_unit.div(by_unit['UNITS'], axis=0).drop(columns=['UNITS'])
    
    by_unit.to_csv(snakemake.output.res_energy_expenses)