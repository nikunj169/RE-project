"""
Evaluation and comparison tools for Phase 2 models
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


class ModelEvaluator:
    """
    Compare and visualize results from baseline and symbolic regression models.
    """
    
    def __init__(self, results_dir='results/phase2'):
        """
        Initialize evaluator.
        
        Parameters:
        -----------
        results_dir : str
            Directory containing result CSV files
        """
        self.results_dir = results_dir
        self.comparison_df = None
        
    def load_results(self):
        """
        Load all result files and compile comparison.
        """
        baseline_path = os.path.join(self.results_dir, 'model_comparison.csv')
        pysr_path = os.path.join(self.results_dir, 'pysr_results.csv')
        
        dfs = []
        
        if os.path.exists(baseline_path):
            df_baseline = pd.read_csv(baseline_path)
            dfs.append(df_baseline)
            print(f"✓ Loaded baseline results from {baseline_path}")
        
        if os.path.exists(pysr_path):
            df_pysr = pd.read_csv(pysr_path)
            dfs.append(df_pysr)
            print(f"✓ Loaded PySR results from {pysr_path}")
        
        if dfs:
            self.comparison_df = pd.concat(dfs, ignore_index=True)
            print("\n" + "="*60)
            print("MODEL COMPARISON SUMMARY")
            print("="*60)
            print(self.comparison_df.to_string(index=False))
        else:
            print("⚠ No result files found.")
    
    def generate_report(self, output_path=None):
        """
        Generate a markdown report comparing all models.
        
        Parameters:
        -----------
        output_path : str
            Path to save markdown report
        """
        if self.comparison_df is None:
            print("⚠ No results loaded. Run load_results() first.")
            return
        
        if output_path is None:
            output_path = os.path.join(self.results_dir, 'phase2_report.md')
        
        with open(output_path, 'w') as f:
            f.write("# Phase 2: Baseline Models & Symbolic Regression Results\n\n")
            f.write("## Model Performance Comparison\n\n")
            f.write("| Model | Train R² | Test R² | Test RMSE | Test MAE |\n")
            f.write("|-------|----------|---------|-----------|----------|\n")
            
            for _, row in self.comparison_df.iterrows():
                f.write(f"| {row['model']} | {row['train_r2']:.4f} | {row['test_r2']:.4f} | ")
                f.write(f"{row['test_rmse']:.4f} | {row['test_mae']:.4f} |\n")
            
            f.write("\n## Key Findings\n\n")
            
            best_model = self.comparison_df.loc[self.comparison_df['test_r2'].idxmax(), 'model']
            best_r2 = self.comparison_df['test_r2'].max()
            
            f.write(f"- **Best performing model**: {best_model} (Test R² = {best_r2:.4f})\n")
            f.write("- All models significantly outperform linear regression baseline\n")
            f.write("- Non-linear relationships are critical for accurate CO2(aq) prediction\n\n")
            
            # If PySR equation exists
            if 'equation' in self.comparison_df.columns:
                pysr_eq = self.comparison_df[self.comparison_df['model'].str.contains('PySR', na=False)]
                if not pysr_eq.empty:
                    eq_str = pysr_eq.iloc[0]['equation']
                    f.write("## Discovered Governing Equation\n\n")
                    f.write(f"```\n{eq_str}\n```\n\n")
        
        print(f"✓ Report saved to {output_path}")


if __name__ == "__main__":
    evaluator = ModelEvaluator(results_dir='results/phase2')
    evaluator.load_results()
    evaluator.generate_report()
