"""E10: sample-matched comparison of TRACE Cant rates with two eMLR(C*) products (GAP_AUDIT_ROUND2).

Products (read-only, data/external/):
  Gruber et al. (2019)  dCant 1994-2007, standard case DCANT_01, ensemble DCANT_01..14
  Mueller et al. (2023) dCant 1994-2004 and 2004-2014, standard case (mask "3", "Standard case"),
                        ensemble = the 11 defined (configuration, mask) cases
For every basin x layer group of E02 (strict QC, same inclusion rules), at the group's samples:
  eMLR rate  = mean of product dCant / period length (nearest horizontal grid point, linear in depth)
  TRACE rate = mean of [Cant(t1) - Cant(t0)] / (t1 - t0) at fixed sample inputs (default TRACE settings)
Estimator-independent: neither rate uses the within-cell regression. Also reports observed beta_DIC and
TRACE beta_Cant from the E02 primary file for context.

Usage: python -m research.experiments.e10_emlr_products
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from research.data_io import ROOT, groups, load_config, load_data
from research.experiments.e02a_run_trace import run_trace

EXT = ROOT / "data" / "external"


class Sampler:
    """Nearest horizontal grid point, linear interpolation in depth; fields held in memory."""

    def __init__(self, lat_grid, lon_grid, dep_grid, lat, lon, depth):
        lat_grid, lon_grid, dep_grid = (np.asarray(x, float) for x in (lat_grid, lon_grid, dep_grid))
        self.i = np.abs(lat[:, None] - lat_grid[None]).argmin(1)
        dl = np.abs(((lon[:, None] - lon_grid[None]) + 180) % 360 - 180)
        self.j = dl.argmin(1)
        k = np.clip(np.searchsorted(dep_grid, depth) - 1, 0, len(dep_grid) - 2)
        self.k0, self.k1 = k, k + 1
        f = (depth - dep_grid[k]) / (dep_grid[k + 1] - dep_grid[k])
        self.f = np.clip(f, 0, 1)

    def __call__(self, field: np.ndarray) -> np.ndarray:  # field dims (depth, lat, lon)
        a = field[self.k0, self.i, self.j]
        b = field[self.k1, self.i, self.j]
        return (1 - self.f) * a + self.f * b


def main() -> None:
    cfg = load_config()
    df = load_data(cfg)
    prim = pd.read_csv(ROOT / "research" / "results" / "e02_layer_trends_primary.csv").set_index(["basin", "layer"])
    g19 = xr.open_dataset(EXT / "gruber2019" / "dcant_emlr_cstar_gruber_94-07_vs1.nc")
    m23 = xr.open_dataset(EXT / "mueller2023" / "dcant_emlr_cstar_mueller_94-04-14.nc")
    g19_fields = {k: g19[f"DCANT_{k:02d}"].transpose("DEPTH", "LATITUDE", "LONGITUDE").values for k in range(1, 15)}
    m_cases = []
    for c in m23.MLR_configuration.values:
        for b in m23.MLR_basin_mask.values:
            sl = m23.dcant.sel(MLR_configuration=c, MLR_basin_mask=b, period="1994 - 2004")
            if np.isfinite(sl.isel(depth=5).values).any():
                m_cases.append((c, b))
    print("Mueller defined cases:", m_cases, flush=True)
    m_fields = {(c, b, per): m23.dcant.sel(MLR_configuration=c, MLR_basin_mask=b, period=per)
                .transpose("depth", "lat", "lon").values
                for c, b in m_cases for per in ("1994 - 2004", "2004 - 2014")}

    rows = []
    for basin, layer, g in groups(df, cfg):
        lat = g.latitude.to_numpy()
        lon = np.where(g.longitude > 180, g.longitude - 360, g.longitude)
        dep = g.depth.to_numpy()
        rec = dict(basin=basin, layer=layer, n=len(g))
        s19 = Sampler(g19.LATITUDE, g19.LONGITUDE, g19.DEPTH, lat, np.mod(lon, 360), dep)
        s23 = Sampler(m23.lat, m23.lon, m23.depth, lat, lon, dep)
        # TRACE fixed-input rates
        tr = {y: np.asarray(run_trace(lon, lat, dep, g.salinity, g.temperature, np.full(len(g), float(y)), cfg)["canth"], float)
              for y in (1994, 2004, 2007, 2014)}
        periods = {"g19_94_07": (1994, 2007), "m23_94_04": (1994, 2004), "m23_04_14": (2004, 2014)}
        for key, (t0, t1) in periods.items():
            rec[f"trace_{key}"] = float(np.nanmean((tr[t1] - tr[t0]) / (t1 - t0)))
        # Gruber 2019
        ens = []
        for k in range(1, 15):
            v = s19(g19_fields[k])
            ens.append(np.nanmean(v) / 13.0)
            if k == 1:
                rec["frac_valid_g19"] = float(np.isfinite(v).mean())
        rec["emlr_g19_94_07"], rec["emlr_g19_min"], rec["emlr_g19_max"] = ens[0], min(ens), max(ens)
        # Mueller 2023
        for per, key in (("1994 - 2004", "m23_94_04"), ("2004 - 2014", "m23_04_14")):
            vals = {}
            for c, b in m_cases:
                v = s23(m_fields[(c, b, per)])
                vals[(c, b)] = np.nanmean(v) / 10.0
                if (c, b) == ("Standard case", "3"):
                    rec[f"frac_valid_{key}"] = float(np.isfinite(v).mean())
            rec[f"emlr_{key}"] = vals[("Standard case", "3")]
            rec[f"emlr_{key}_min"], rec[f"emlr_{key}_max"] = min(vals.values()), max(vals.values())
        p = prim.loc[(basin, layer)]
        rec.update(beta_dic=p.beta_tco2, beta_cant=p.beta_canth, R_obs=p.R)
        for key in periods:
            src = "g19_94_07" if key.startswith("g19") else key
            rec[f"ratio_emlr_trace_{key}"] = rec[f"emlr_{src}"] / rec[f"trace_{key}"]
        rows.append(rec)
        print(basin, layer, {k: round(v, 3) for k, v in rec.items() if k.startswith("ratio") or k == "R_obs"}, flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "research" / "results" / "e10_emlr_products.csv", index=False)


if __name__ == "__main__":
    main()
