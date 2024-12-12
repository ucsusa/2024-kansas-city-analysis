import pandas as pd
import numpy as np
import pypsa
import matplotlib.pyplot as plt


if __name__ == "__main__":
    solver = snakemake.config['solver']

    n = pypsa.Network(snakemake.input.elec_network)

    try:
        n.optimize(solver_name=solver)
    except:
        print(f'No {solver} solver found. Reverting to HiGHS')
        n.optimize(solver_name='highs')

    n.export_to_netcdf(snakemake.output.residential_model_solved)