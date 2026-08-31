import numpy as np
import pandas as pd
from scipy.io import loadmat
from sklearn.metrics import mean_squared_error
from datetime import datetime
from canyonbpy import canyonb

print("Loading MAT file...")

data = loadmat("data/raw/GLODAPv2.2023_Merged_Master_File.mat")

lat = data["G2latitude"].flatten()
lon = data["G2longitude"].flatten()
pres = data["G2pressure"].flatten()

temp = data["G2temperature"].flatten()
sal = data["G2salinity"].flatten()
oxy = data["G2oxygen"].flatten()

dic_true = data["G2tco2"].flatten()

print("Building dataframe...")

df = pd.DataFrame({
    "lat": lat,
    "lon": lon,
    "pres": pres,
    "temp": temp,
    "sal": sal,
    "oxy": oxy,
    "dic_true": dic_true
})

df = df.dropna()

print("Rows after cleaning:", len(df))

print("Running CANYON-B predictions...")

results = canyonb(
    gtime=[datetime(2020,1,1)] * len(df),
    lat=df["lat"].values,
    lon=df["lon"].values,
    pres=df["pres"].values,
    temp=df["temp"].values,
    psal=df["sal"].values,
    doxy=df["oxy"].values
)

df["canyon_dic"] = results["CT"]

print("Computing RMSE...")

rmse = np.sqrt(mean_squared_error(df["dic_true"], df["canyon_dic"]))

print("\nCANYON-B RMSE:", rmse)

df.to_csv("glodap_canyonb_results.csv", index=False)

print("Saved predictions")