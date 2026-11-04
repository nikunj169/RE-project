"""Build the research dataset from the raw GLODAPv2.2023 merged master file.

The raw file is read-only and never modified. Output: research/data/glodap_research.parquet
(or .csv.gz fallback) with all variables needed by the research program plus QC tiers.

QC tiers (column ``qc_tier``):
  strict  : G2tco2f==2 & G2tco2qc==1 & salinity/oxygen WOCE flag 2 & salinity/oxygen qc==1
  relaxed : the repository's original range-only QC (for ablation / continuity)

Usage:
    python research/preprocessing/build_dataset.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "GLODAPv2.2023_Merged_Master_File.mat"
OUT_DIR = ROOT / "research" / "data"
EXPECTED_SHA256 = "002881fa71923d15c1bd5aca3b99bc220ad5d1c77a90f130ba1ff0d3910d8766"

VARS = {
    "cruise": "G2cruise", "station": "G2station", "cast": "G2cast", "region": "G2region",
    "year": "G2year", "month": "G2month", "day": "G2day",
    "latitude": "G2latitude", "longitude": "G2longitude", "bottomdepth": "G2bottomdepth",
    "pressure": "G2pressure", "depth": "G2depth", "temperature": "G2temperature", "theta": "G2theta",
    "salinity": "G2salinity", "salinity_f": "G2salinityf", "salinity_qc": "G2salinityqc",
    "gamma": "G2gamma", "sigma0": "G2sigma0",
    "oxygen": "G2oxygen", "oxygen_f": "G2oxygenf", "oxygen_qc": "G2oxygenqc", "aou": "G2aou",
    "nitrate": "G2nitrate", "nitrate_f": "G2nitratef", "phosphate": "G2phosphate", "phosphate_f": "G2phosphatef",
    "silicate": "G2silicate", "silicate_f": "G2silicatef",
    "tco2": "G2tco2", "tco2_f": "G2tco2f", "tco2_qc": "G2tco2qc",
    "talk": "G2talk", "talk_f": "G2talkf", "talk_qc": "G2talkqc",
    "cfc12": "G2cfc12", "pcfc12": "G2pcfc12", "cfc12_f": "G2cfc12f",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def assign_basin(lat: np.ndarray, region: np.ndarray) -> np.ndarray:
    """Same definition as src/basins.py, plus explicit Arctic label instead of 'Unclassified'."""
    out = np.full(lat.shape, "Unclassified", dtype=object)
    out[lat < -35] = "Southern Ocean"
    for code, name in ((1, "Atlantic"), (16, "Indian"), (8, "Pacific")):
        out[(region == code) & (lat >= -35)] = name
    out[(region == 4) & (lat >= -35)] = "Arctic"
    return out


def main() -> None:
    digest = sha256(RAW)
    if digest != EXPECTED_SHA256:
        raise RuntimeError(f"Raw file checksum mismatch: {digest}")
    mat = sio.loadmat(RAW, squeeze_me=True)
    df = pd.DataFrame({k: np.asarray(mat[v], dtype=float) for k, v in VARS.items()})
    expocodes = np.asarray(mat["expocode"]).astype(str)
    expono = np.asarray(mat["expocodeno"], dtype=float)
    df["expocode"] = pd.Series(df["cruise"]).map(dict(zip(expono, expocodes)))
    n_raw = len(df)

    core = df[["salinity", "temperature", "aou", "tco2"]].notna().all(axis=1)
    rng = ((df.salinity > 25) & (df.salinity < 42) & (df.temperature > -2.5) & (df.temperature < 35)
           & (df.tco2 > 1700) & (df.tco2 < 2600) & (df.aou > -50))
    relaxed = core & rng
    strict = (relaxed & (df.tco2_f == 2) & (df.tco2_qc == 1)
              & (df.salinity_f == 2) & (df.salinity_qc == 1)
              & (df.oxygen_f == 2) & (df.oxygen_qc == 1) & df.theta.notna() & df.pressure.notna())

    out = df[relaxed].copy()
    out["qc_strict"] = strict[relaxed].to_numpy()
    out["basin"] = assign_basin(out.latitude.to_numpy(), out.region.to_numpy())
    out["decimal_year"] = out.year + (out.month.fillna(7) - 0.5) / 12.0
    out["profile_id"] = out.cruise.astype(int).astype(str) + "_" + out.station.astype(int).astype(str)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / "glodap_research.csv.gz"
    out.to_csv(path, index=False)
    stats = {
        "raw_sha256": digest,
        "n_raw": int(n_raw),
        "n_relaxed_repo_qc": int(relaxed.sum()),
        "n_strict": int(strict.sum()),
        "dropped_by_strict": {
            "tco2_calculated_or_flagged": int((relaxed & (df.tco2_f != 2)).sum()),
            "tco2_no_secondary_qc": int((relaxed & (df.tco2_qc != 1)).sum()),
            "oxygen_or_salinity_flag_or_qc": int((relaxed & ~((df.salinity_f == 2) & (df.salinity_qc == 1)
                                                              & (df.oxygen_f == 2) & (df.oxygen_qc == 1))).sum()),
        },
        "output": str(path.relative_to(ROOT)),
    }
    (OUT_DIR / "build_stats.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
