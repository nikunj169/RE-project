"""
Main script to run complete Phase 2 pipeline:
1. Train baseline models (RF, XGBoost)
2. Run symbolic regression (PySR)
3. Generate comparison report
"""

import sys
import os
import argparse

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from phase2.baseline_models import BaselineModelTrainer
from phase2.symbolic_regression import SymbolicRegressionDiscovery
from phase2.evaluation import ModelEvaluator


def run_full_phase2(data_path, 
                    skip_xgboost=False, 
                    skip_pysr=False,
                    pysr_sample_frac=0.5,
                    pysr_iterations=40):
    """
    Run complete Phase 2 analysis pipeline.
    
    Parameters:
    -----------
    data_path : str
        Path to processed CSV from Phase 1
    skip_xgboost : bool
        Skip XGBoost training (if not installed)
    skip_pysr : bool
        Skip symbolic regression (if PySR not installed or for quick testing)
    pysr_sample_frac : float
        Fraction of data to use for PySR (for computational efficiency)
    pysr_iterations : int
        Number of PySR iterations
    """
    
    print("\n" + "="*70)
    print(" PHASE 2: BASELINE MODELS & SYMBOLIC REGRESSION")
    print("="*70)
    print(f"\nData source: {data_path}\n")
    
    # =====================================================================
    # STEP 1: Train Baseline Models
    # =====================================================================
    print("\n" + "="*70)
    print(" STEP 1: TRAINING BASELINE MODELS")
    print("="*70 + "\n")
    
    trainer = BaselineModelTrainer(
        data_path=data_path,
        test_size=0.2,
        random_state=42
    )
    
    trainer.load_and_prepare_data()
    
    # Train Random Forest
    trainer.train_random_forest(n_estimators=200, max_depth=None)
    
    # Train XGBoost (if available and not skipped)
    if not skip_xgboost:
        trainer.train_xgboost(n_estimators=200, max_depth=6)
    
    # Compute SHAP values for Random Forest
    trainer.compute_shap_values('RandomForest', sample_size=500)
    
    # Save models and results
    trainer.save_models(save_dir='models/phase2')
    trainer.save_results(save_dir='results/phase2')
    
    # =====================================================================
    # STEP 2: Symbolic Regression (Optional)
    # =====================================================================
    if not skip_pysr:
        print("\n" + "="*70)
        print(" STEP 2: SYMBOLIC REGRESSION (PySR)")
        print("="*70 + "\n")
        
        sr = SymbolicRegressionDiscovery(
            data_path=data_path,
            test_size=0.2,
            random_state=42
        )
        
        sr.load_and_prepare_data(sample_frac=pysr_sample_frac)
        
        sr.run_symbolic_regression(
            populations=15,
            population_size=50,
            niterations=pysr_iterations,
            maxsize=20,
            binary_operators=["+", "-", "*", "/"],
            unary_operators=["exp", "log"]
        )
        
        sr.save_model(save_dir='models/phase2')
        sr.save_results(save_dir='results/phase2')
    else:
        print("\n⚠ Skipping symbolic regression (--skip-pysr flag)")
    
    # =====================================================================
    # STEP 3: Generate Comparison Report
    # =====================================================================
    print("\n" + "="*70)
    print(" STEP 3: GENERATING COMPARISON REPORT")
    print("="*70 + "\n")
    
    evaluator = ModelEvaluator(results_dir='results/phase2')
    evaluator.load_results()
    evaluator.generate_report(output_path='results/phase2/phase2_report.md')
    
    # =====================================================================
    # COMPLETION
    # =====================================================================
    print("\n" + "="*70)
    print(" PHASE 2 COMPLETE!")
    print("="*70)
    print("\n✓ All models trained and evaluated")
    print("✓ Results saved to: results/phase2/")
    print("✓ Models saved to:  models/phase2/")
    print("\nNext steps:")
    print("  1. Review results/phase2/phase2_report.md")
    print("  2. Examine feature importances in results/phase2/")
    print("  3. Analyze discovered equation in models/phase2/discovered_equation.txt")
    print("  4. Proceed to Phase 3 for physical interpretation\n")


def main():
    """
    Command-line interface for Phase 2.
    """
    parser = argparse.ArgumentParser(
        description="Phase 2: Baseline Models & Symbolic Regression"
    )
    
    parser.add_argument(
        '--data',
        type=str,
        default='data/processed/southern_ocean_training.csv',
        help='Path to processed CSV from Phase 1'
    )
    
    parser.add_argument(
        '--skip-xgboost',
        action='store_true',
        help='Skip XGBoost training'
    )
    
    parser.add_argument(
        '--skip-pysr',
        action='store_true',
        help='Skip symbolic regression (faster for testing)'
    )
    
    parser.add_argument(
        '--pysr-sample',
        type=float,
        default=0.5,
        help='Fraction of data for PySR (0.1-1.0, default=0.5)'
    )
    
    parser.add_argument(
        '--pysr-iterations',
        type=int,
        default=40,
        help='Number of PySR iterations (default=40)'
    )
    
    args = parser.parse_args()
    
    # Validate data file exists
    if not os.path.exists(args.data):
        print(f"❌ ERROR: Data file not found: {args.data}")
        print("\nMake sure you've run Phase 1 preprocessing first!")
        sys.exit(1)
    
    # Run Phase 2
    run_full_phase2(
        data_path=args.data,
        skip_xgboost=args.skip_xgboost,
        skip_pysr=args.skip_pysr,
        pysr_sample_frac=args.pysr_sample,
        pysr_iterations=args.pysr_iterations
    )


if __name__ == "__main__":
    main()
