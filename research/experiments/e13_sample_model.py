"""E13/E15 step 1: sample a CMIP6 Earth system model at every E02 sample (GAP_AUDIT_ROUND3 §1, §8).

Annual Oyr dissic, dissicnat, o2 and monthly Omon thetao, so (averaged to calendar years in memory)
are read lazily from the public Pangeo CMIP6 Zarr store (historical 1972-2014, ssp585 2015-2021,
member r1i2p1f1, regular 1-degree grid "gr"). Each sample takes its own year; horizontally the
nearest model column that is wet at the sample's depth, linear in depth.

Models (MODELS): MRI-ESM2-0 (E13; regular grid) and NorESM2-LM (E15; curvilinear grid with 2-D
lat/lon, historical only, so samples after 2014 stay NaN).

Output: research/data/model_<name>.csv.gz (row_id; model theta, salinity, temperature, aou,
oxygen, tco2, tco2nat, canth_true; all per kg via GSW in-situ density).

Usage: python -m research.experiments.e13_sample_model [--model mri_esm2_0|noresm2_lm]
"""

from __future__ import annotations

import argparse
import time
from concurrent.futures import ThreadPoolExecutor

import gsw
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree

from research.data_io import ROOT, groups, load_config, load_data

CATALOG = "https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv"
MODELS = {"mri_esm2_0": ("MRI-ESM2-0", "r1i2p1f1", ("historical", "ssp585")),
          "noresm2_lm": ("NorESM2-LM", "r2i1p1f1", ("historical",))}
VARS = {"dissic": "Oyr", "dissicnat": "Oyr", "o2": "Oyr", "thetao": "Omon", "so": "Omon"}


def open_stores(model, member, experiments):
    local = ROOT / "research" / "data" / "_cache_cmip6_catalog.csv"  # ~100 MB; download once
    if not local.exists():
        pd.read_csv(CATALOG).to_csv(local, index=False)
    cat = pd.read_csv(local)
    cat = cat[(cat.source_id == model) & (cat.member_id == member) & (cat.grid_label == "gr")]
    stores = {}
    for exp in experiments:
        for v, tab in VARS.items():
            r = cat[(cat.experiment_id == exp) & (cat.variable_id == v) & (cat.table_id == tab)]
            assert len(r) == 1, (exp, v, len(r))
            r = r.iloc[0]
            print(exp, v, tab, r.version, r.zstore, flush=True)
            stores[exp, v] = xr.open_zarr(r.zstore, consolidated=True, storage_options={"token": "anon"})[v]
    return stores


def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


class WetSampler:
    """Nearest model column wet at both bracketing levels; linear interpolation in depth.

    lat_g/lon_g are 1-D (regular grid) or 2-D (curvilinear, same shape as the horizontal field)."""

    def __init__(self, lat_g, lon_g, dep_g, wet, lat, lon, depth):
        kmax = int(np.nonzero(wet.any(axis=(1, 2)))[0].max())  # deepest level with any ocean
        k = np.clip(np.searchsorted(dep_g, depth) - 1, 0, kmax - 1)
        self.k0, self.k1 = k, k + 1
        self.f = np.clip((depth - dep_g[k]) / (dep_g[k + 1] - dep_g[k]), 0, 1)
        LA, LO = np.meshgrid(lat_g, lon_g, indexing="ij") if np.ndim(lat_g) == 1 else (lat_g, lon_g)
        self.i = np.empty(len(lat), int)
        self.j = np.empty(len(lat), int)
        pts = xyz(lat, lon)
        for kk in np.unique(self.k1):
            sel = self.k1 == kk
            ii, jj = np.nonzero(wet[kk])
            _, nn = cKDTree(xyz(LA[ii, jj], LO[ii, jj])).query(pts[sel])
            self.i[sel], self.j[sel] = ii[nn], jj[nn]

    def __call__(self, field, rows):  # field (lev, lat, lon)
        a = field[self.k0[rows], self.i[rows], self.j[rows]]
        b = field[self.k1[rows], self.i[rows], self.j[rows]]
        return (1 - self.f[rows]) * a + self.f[rows] * b


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="mri_esm2_0", choices=list(MODELS))
    args = ap.parse_args()
    model, member, experiments = MODELS[args.model]
    out_file = ROOT / "research" / "data" / f"model_{args.model}.csv.gz"
    cfg = load_config()
    df = load_data(cfg)
    u = pd.concat([g for _, _, g in groups(df, cfg)])
    print("samples", len(u), flush=True)
    stores = open_stores(model, member, experiments)
    ref = stores["historical", "thetao"]
    wet = np.isfinite(ref.isel(time=0).values)
    lon = np.mod(u.longitude.to_numpy(float), 360)
    lat_g, lon_g = (ref.lat.values, ref.lon.values) if "lat" in ref.coords else (ref.latitude.values, ref.longitude.values)
    samp = WetSampler(lat_g, lon_g, ref.lev.values, wet,
                      u.latitude.to_numpy(float), lon, u.depth.to_numpy(float))
    year = u.year.to_numpy(int)
    out = {v: np.full(len(u), np.nan) for v in VARS}
    cache = ROOT / "research" / "data" / f"_cache_model_{args.model}"  # per-year checkpoints (resume after a hung read)
    cache.mkdir(exist_ok=True)

    def annual(exp, v, y):
        da = stores[exp, v]
        idx = np.nonzero(da.time.dt.year.values == y)[0]
        assert len(idx) == (12 if VARS[v] == "Omon" else 1), (v, y, len(idx))
        return v, da.isel(time=idx).values.mean(0)

    t0 = time.time()
    with ThreadPoolExecutor(5) as pool:
        for y in range(year.min(), year.max() + 1):
            rows = np.nonzero(year == y)[0]
            if not len(rows):
                continue
            exp = "historical" if y <= 2014 else "ssp585"
            if exp not in experiments:
                continue
            ck = cache / f"{y}.npz"
            if ck.exists():
                z = np.load(ck)
                assert np.array_equal(z["rows"], rows)
                for v in VARS:
                    out[v][rows] = z[v]
                print(f"{y}: {len(rows)} samples (cached)", flush=True)
                continue
            for v, field in pool.map(lambda v: annual(exp, v, y), VARS):
                out[v][rows] = samp(field, rows)
            np.savez(ck, rows=rows, **{v: out[v][rows] for v in VARS})
            print(f"{y}: {len(rows)} samples, {time.time() - t0:.0f}s", flush=True)

    p = u.pressure.to_numpy(float)
    lat = u.latitude.to_numpy(float)
    SA = gsw.SA_from_SP(out["so"], p, lon, lat)
    CT = gsw.CT_from_pt(SA, out["thetao"])
    rho = gsw.rho(SA, CT, p)
    k = 1e6 / rho  # mol m-3 -> umol kg-1
    o2 = out["o2"] * k
    res = pd.DataFrame({
        "row_id": u.index.to_numpy(),
        "theta": out["thetao"], "salinity": out["so"], "temperature": gsw.t_from_CT(SA, CT, p),
        "oxygen": o2, "aou": gsw.O2sol_SP_pt(out["so"], out["thetao"]) - o2,
        "tco2": out["dissic"] * k, "tco2nat": out["dissicnat"] * k,
    })
    res["canth_true"] = res.tco2 - res.tco2nat
    print(res.describe().T.round(2).to_string(), flush=True)
    print("non-finite rows:", int((~np.isfinite(res.drop(columns="row_id"))).any(axis=1).sum()))
    res.to_csv(out_file, index=False, float_format="%.5f")
    print("wrote", out_file)


if __name__ == "__main__":
    main()
