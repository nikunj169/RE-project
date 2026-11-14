"""E06: drift-immune differential trends (target layer relative to deep reference, cruise FE).

Pre-registered in research_docs/GAP_AUDIT_ROUND2.md. For each basin and each target layer with
upper edge <= 27.8, the target is stacked with the reference layers (lower edge >= 27.8) using the
same inclusion rules as E02, restricted to cruises sampling both. DIC and TRACE Cant are fitted
with the identical two-way FE estimator; Ddelta = dbeta(DIC) - dbeta(Cant).

Also writes the naive comparison delta_target - delta_ref from single-layer E02 fits on the same rows.

Usage:
  python -m research.experiments.e06_differential_trends [--boot 300]
  python -m research.experiments.e06_differential_trends --synthetic [--reps 40]
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from research.data_io import ROOT, load_config, load_data
from research.models.differential_trend import _codes, design, exact_twoway, fit_differential
from research.models.within_cell_trend import fit_trend, select_cells

SEED = 20261114
REF_MIN = 27.8


def lower(layer: str) -> float:
    return float(layer.strip("(]").split(",")[0])


def upper(layer: str) -> float:
    return float(layer.strip("(]").split(",")[1])


def basin_frames(cfg: dict, cant_file: str = "trace_canth.csv.gz", df: pd.DataFrame | None = None,
                 ref: tuple[float, float] = (REF_MIN, np.inf), target_max: float = REF_MIN):
    """Yield (basin, target_layer, stacked frame, target mask) following E02 inclusion rules.

    If df is given it is used as is (it must already carry a canth column). Reference layers are
    those with ref[0] <= lower edge < ref[1]; targets are layers with upper edge <= target_max."""
    if df is None:
        df = load_data(cfg)
        tr = pd.read_csv(ROOT / "research" / "data" / cant_file).set_index("row_id")
        df = df.join(tr[["canth"]], how="inner")
    df = df[np.isfinite(df.canth)]
    for basin, b in df.groupby("basin"):
        kept = {}
        for layer, g in b.groupby("layer"):
            g = select_cells(g, cfg["min_span_years"])
            if len(g) >= cfg["min_samples"] and g["cell"].nunique() >= cfg["min_cells"]:
                kept[layer] = g
        ref_df = pd.concat([g for l, g in kept.items() if ref[0] <= lower(l) < ref[1]])
        for layer, g in sorted(kept.items(), key=lambda kv: lower(kv[0])):
            if upper(layer) > target_max:
                continue
            both = np.intersect1d(g.cruise.unique(), ref_df.cruise.unique())
            st = pd.concat([g[g.cruise.isin(both)], ref_df[ref_df.cruise.isin(both)]])
            mask = np.r_[np.ones(g.cruise.isin(both).sum(), bool), np.zeros(ref_df.cruise.isin(both).sum(), bool)]
            yield basin, layer, st, mask


def run_observed(cfg: dict, n_boot: int) -> None:
    rows = []
    t0 = time.time()
    for basin, layer, st, mask in basin_frames(cfg):
        fit = fit_differential(st, mask, ["tco2", "canth"], n_boot=n_boot, seed=SEED)
        bd, bc = fit["boot"]["tco2"], fit["boot"]["canth"]
        dd = bd - bc
        # naive: single-layer E02 estimator on the same rows, target minus pooled reference
        ft = fit_trend(st[mask], ["tco2", "canth"])
        fr = fit_trend(st[~mask], ["tco2", "canth"])
        naive = (ft.beta["tco2"] - ft.beta["canth"]) - (fr.beta["tco2"] - fr.beta["canth"])
        q = lambda x: np.nanquantile(x, [0.025, 0.975])
        rec = dict(basin=basin, layer=layer, n=fit["n"], n_target=fit["n_target"], n_cruises=fit["n_cruises"],
                   dbeta_dic=fit["point"]["tco2"], dbeta_cant=fit["point"]["canth"],
                   ddelta=fit["point"]["tco2"] - fit["point"]["canth"],
                   ddelta_lo=q(dd)[0], ddelta_hi=q(dd)[1],
                   dbeta_dic_lo=q(bd)[0], dbeta_dic_hi=q(bd)[1],
                   delta_target_same_rows=ft.beta["tco2"] - ft.beta["canth"],
                   delta_ref_same_rows=fr.beta["tco2"] - fr.beta["canth"], naive_ddelta=naive,
                   n_boot_ok=int(np.isfinite(dd).sum()))
        rows.append(rec)
        print(f"{basin:15s} {layer:13s} n={rec['n']:6d} cru={rec['n_cruises']:3d} "
              f"Ddelta={rec['ddelta']:+.3f} [{rec['ddelta_lo']:+.3f},{rec['ddelta_hi']:+.3f}] "
              f"naive={naive:+.3f} delta_ref={rec['delta_ref_same_rows']:+.3f} ({time.time()-t0:.0f}s)", flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "research" / "results" / "e06_differential_trends.csv", index=False)


def run_synthetic(cfg: dict, reps: int, n_boot: int) -> None:
    """Bias of E02 (single-layer) vs E06 (differential) under time-correlated cruise drift."""
    targets = {("Pacific", "(26.0, 26.5]"), ("Atlantic", "(26.5, 27.0]"), ("Indian", "(26.5, 27.0]")}
    rng = np.random.default_rng(SEED)
    rows = []
    for basin, layer, st, mask in basin_frames(cfg):
        if (basin, layer) not in targets:
            continue
        st = st.copy()
        # static null component: two-way FE + conditioning without time
        X = design(st, mask, ("theta", "salinity", "aou"))[:, :-1]
        a = _codes(st["cell"].astype(str) + "|" + np.where(mask, "T", "R"))
        cru = _codes(st["cruise"])
        y = st["tco2"].to_numpy(float)
        _, resid = exact_twoway(X, y[:, None], a, cru, np.ones(len(y)), return_resid=True)
        resid = resid[:, 0]
        fitted = y - resid
        t = st["decimal_year"].to_numpy() - 2000
        tc = st.groupby("cruise")["decimal_year"].transform("mean").to_numpy() - 2000
        b_t, b_r = 0.8, 0.1
        for d in (0.0, 0.05, 0.1):
            e02, e06, cover = [], [], []
            for r in range(reps):
                off = (d * tc + rng.normal(0, 3.0, cru.max() + 1)[cru])
                st["y_syn"] = fitted + np.where(mask, b_t, b_r) * t + off + rng.permutation(resid)
                e02.append(fit_trend(st[mask], ["y_syn"]).beta["y_syn"] - fit_trend(st[~mask], ["y_syn"]).beta["y_syn"])
                f = fit_differential(st, mask, ["y_syn"], n_boot=n_boot, seed=int(rng.integers(1 << 31)))
                e06.append(f["point"]["y_syn"])
                if n_boot:
                    lo, hi = np.nanquantile(f["boot"]["y_syn"], [0.025, 0.975])
                    cover.append(lo <= b_t - b_r <= hi)
            single = [fit_trend(st[mask].assign(y_syn=fitted[mask] + b_t * t[mask] + d * tc[mask]), ["y_syn"]).beta["y_syn"]]
            rows.append(dict(basin=basin, layer=layer, drift=d, true_ddelta=b_t - b_r,
                             e02_target_bias_noiseless=single[0] - b_t,
                             e02_diff_bias=float(np.mean(e02) - (b_t - b_r)),
                             e06_bias=float(np.mean(e06) - (b_t - b_r)), e06_sd=float(np.std(e06)),
                             e06_coverage=float(np.mean(cover)) if cover else np.nan, reps=reps))
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "research" / "results" / "e06_synthetic_drift.csv", index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=300)
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--reps", type=int, default=40)
    args = ap.parse_args()
    cfg = load_config()
    if args.synthetic:
        run_synthetic(cfg, args.reps, args.boot)
    else:
        run_observed(cfg, args.boot)


if __name__ == "__main__":
    main()
