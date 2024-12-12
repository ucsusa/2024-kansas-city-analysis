import geopandas as gpd
import os
import time
import zipfile
import requests


if __name__ == "__main__":

    if not os.path.exists("..\\..\\cjest-data\\"):
        print('Downloading file from internet')
        r = requests.get(cjest_url)
        z = zipfile.ZipFile(io.BytesIO(r.content))
        z.extractall("..\\..\\cjest-data\\")

        start = time.perf_counter()
        cjest_df = gpd.read_file("..\\..\\cjest-data\\usa.zip")
        end = time.perf_counter()
        print(f"It took {(end-start)/60:3f} minutes to load the CJEST data")
    else:
        print("Loading from previously saved file...")
        start = time.perf_counter()
        cjest_df = gpd.read_file("..\\..\\cjest-data\\usa.zip")
        end = time.perf_counter()
        print(f"It took {(end-start)/60:3f} minutes to load the CJEST data")