"""Cross-check core manuscript claims against generated result files."""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "FINAL_MANUSCRIPT_NUMERIC_AUDIT.md"


def main():
    processed = ROOT / "results" / "processed"
    raw = ROOT / "results" / "raw" / "symbolic_regression"
    sr = pd.read_csv(raw / "symbolic_candidates_all.csv")
    recurrence = pd.read_csv(ROOT / "results" / "tables" / "equation_family_recurrence.csv")
    hierarchy = pd.read_csv(processed / "final_model_hierarchy.csv")
    blocked = pd.read_csv(processed / "spatial_cruise_validation_summary.csv")
    watermass = pd.read_csv(processed / "watermass_atlantic.csv")
    abyss = pd.read_csv(processed / "southern_ocean_abyssal.csv").iloc[0]

    lines = ["# Final Manuscript Numeric Audit", "", "All checks below are computed from generated result files.", ""]
    lines.append(f"- Symbolic candidates: {len(sr)}")
    lines.append(f"- Rank-1 candidates: {(sr['family'] == 'rank-1 quadratic').sum()}")
    lines.append(f"- Rank-1 cells: {recurrence['rank1_present'].sum()}/{len(recurrence)}")
    lines.append(f"- Spatial/cruise validation rows: {len(blocked)}")
    lines.append("")

    lines.append("## Rank-1 temporal and external RMSE")
    r1 = hierarchy[hierarchy.model == "rank1"]
    lines.append(r1[["basin", "validation_RMSE", "external_RMSE", "parameter_count", "expression_complexity"]].to_markdown(index=False))
    lines.append("")

    lines.append("## Rank-2 deltas")
    rank2 = hierarchy[hierarchy.model == "rank2"].merge(
        r1[["basin", "validation_RMSE", "external_RMSE"]], on="basin", suffixes=("_rank2", "_rank1")
    )
    rank2["validation_delta"] = rank2.validation_RMSE_rank2 - rank2.validation_RMSE_rank1
    rank2["external_delta"] = rank2.external_RMSE_rank2 - rank2.external_RMSE_rank1
    lines.append(rank2[["basin", "validation_delta", "external_delta"]].to_markdown(index=False))
    lines.append("")

    lines.append("## Blocked validation")
    lines.append(blocked[["basin", "block_type", "model", "n_test", "n_groups_total", "RMSE_mean", "RMSE_sd", "R2_mean", "Bias_mean"]].to_markdown(index=False))
    lines.append("")

    lines.append("## AAIW and abyssal failure")
    aaiw = watermass[watermass.water_mass == "AAIW"]
    lines.append(aaiw[["Model", "RMSE", "Bias", "R2", "n"]].to_markdown(index=False))
    lines.append(f"\nSouthern abyss global R2={abyss.global_r2:.6f}, local R2={abyss.local_r2:.6f}, local RMSE={abyss.local_rmse:.6f}.")
    lines.append("")

    lines.append("## Consistency checks")
    checks = {
        "candidate_count_599": len(sr) == 599,
        "rank1_count_86": int((sr["family"] == "rank-1 quadratic").sum()) == 86,
        "rank1_cells_24_of_48": int(recurrence["rank1_present"].sum()) == 24 and len(recurrence) == 48,
        "blocked_rows_pre2018": bool((blocked["n_test"] > 0).all()),
        "cruise_and_spatial_present": set(blocked.block_type) == {"spatial", "cruise"},
        "abyss_global_negative": abyss.global_r2 < 0,
        "abyss_local_positive": abyss.local_r2 > 0.7,
    }
    for key, value in checks.items():
        lines.append(f"- {key}: {'PASS' if value else 'FAIL'}")

    OUT.write_text("\n".join(lines) + "\n")
    failed = [key for key, value in checks.items() if not value]
    if failed:
        raise SystemExit(f"Numeric audit failed: {failed}")
    print(OUT)


if __name__ == "__main__":
    main()
