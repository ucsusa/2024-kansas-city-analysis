import pandas as pd
import numpy as np
from tqdm import tqdm

outage_files = {2018:"https://figshare.com/ndownloader/files/42547879",
                2019:"https://figshare.com/ndownloader/files/42547885",
                2020:"https://figshare.com/ndownloader/files/42547894",
                2021:"https://figshare.com/ndownloader/files/42547891",
                2022:"https://figshare.com/ndownloader/files/42547897",
                2023:"https://figshare.com/ndownloader/files/44574907"}

if __name__ == "__main__":

    frames = []
    pbar = tqdm(outage_files.items())
    for year, url in pbar:
        pbar.set_description(desc=f"{year}")
        df = pd.read_csv(url, parse_dates=True, index_col='run_start_time')
        frames.append(df)

    outages = pd.concat(frames, axis=0)


    outages = outages.loc[outages['county'] == snakemake.config['county'].capitalize()]

    outages.to_csv(snakemake.output.outages)