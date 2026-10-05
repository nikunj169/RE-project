"""
Master pipeline script.

Runs the complete analysis pipeline in sequence:
  1. Preprocessing
  2. Baseline reproduction
  3. Model family comparison
  4. PySR Pareto search (requires Julia + PySR)
  5. Equivalence analysis
  6. Depth/water-mass/Southern Ocean analysis
  7. Figure generation
  8. Table generation

Usage:
    python scripts/run_all.py [--skip-pysr]
"""

import sys
import os
import logging
import argparse
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("logs/pipeline.log"),
    ])
logger = logging.getLogger(__name__)


def run_step(name, script_path, skip=False):
    """Run a pipeline step."""
    if skip:
        logger.info(f"SKIPPING: {name}")
        return

    logger.info(f"\n{'=' * 70}")
    logger.info(f"RUNNING: {name}")
    logger.info(f"{'=' * 70}")

    start = time.time()
    exit_code = os.system(f"python {script_path}")
    elapsed = time.time() - start

    if exit_code != 0:
        logger.error(f"FAILED: {name} (exit code {exit_code}, {elapsed:.1f}s)")
    else:
        logger.info(f"COMPLETED: {name} ({elapsed:.1f}s)")


def main():
    parser = argparse.ArgumentParser(description="Run the full analysis pipeline")
    parser.add_argument("--skip-pysr", action="store_true",
                        help="Skip PySR steps (requires Julia)")
    parser.add_argument("--start-from", type=int, default=1,
                        help="Start from step N")
    args = parser.parse_args()

    os.makedirs("logs", exist_ok=True)

    script_dir = os.path.dirname(os.path.abspath(__file__))

    steps = [
        ("Step 01: Preprocessing", "01_preprocess.py", False),
        ("Step 02: Baseline Reproduction", "02_baseline_reproduction.py", False),
        ("Step 03: Model Families", "03_model_families.py", False),
        ("Step 04: PySR Pareto Search", "04_pysr_pareto_search.py", args.skip_pysr),
        ("Step 05: Equivalence Analysis", "05_equivalence_analysis.py", False),
        ("Step 06: Depth/WM/SO Analysis", "06_depth_watermass_southern_ocean.py", False),
        ("Step 07: Generate Figures", "07_generate_figures.py", False),
        ("Step 08: Generate Tables", "08_generate_tables.py", False),
    ]

    logger.info("=" * 70)
    logger.info("FULL PIPELINE EXECUTION")
    logger.info("=" * 70)

    for i, (name, script, skip) in enumerate(steps, 1):
        if i < args.start_from:
            logger.info(f"Skipping {name} (start_from={args.start_from})")
            continue
        run_step(name, os.path.join(script_dir, script), skip=skip)

    logger.info("\n" + "=" * 70)
    logger.info("PIPELINE COMPLETE")
    logger.info("=" * 70)
    logger.info("Results in: results/processed/")
    logger.info("Figures in: results/figures/")
    logger.info("Tables in:  results/tables/")


if __name__ == "__main__":
    main()
