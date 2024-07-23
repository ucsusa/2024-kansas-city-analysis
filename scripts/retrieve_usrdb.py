import pandas as pd


if __name__ == "__main__":
    
    URL = "https://apps.openei.org/USURDB/download/usurdb.csv.gz"
    usrdb = pd.read_csv(URL, low_memory=False, parse_dates=True)