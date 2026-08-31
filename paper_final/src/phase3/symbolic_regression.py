"""
Phase 3: Symbolic Regression - Discover Governing Equation
Uses PySR to find interpretable mathematical relationship:
CO₂(aq) = f(Salinity, Temperature, AOU, Stratification, C_atm)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.preprocessing import StandardScaler
import os
import warnings
warnings.filterwarnings('ignore')

try:
    from pysr import PySRRegressor
    PYSR_AVAILABLE = True
except ImportError:
    PYSR_AVAILABLE = False
    print("⚠ PySR not installed. Run: pip install pysr && python -m pysr install")

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300


class SymbolicRegressionPhase3:
    """
    Discover interpretable governing equations for CO₂(aq).
    """
    
    def __init__(self, 
                 data_path='data/processed/southern_ocean_training.csv',
                 output_dir='results/phase3'):
        """
        Initialize Phase 3 analysis.
        """
        self.data_path = data_path
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        self.df = None
        self.pysr_models = {}
        self.results = {}
        
    def load_data(self):
        """Load and prepare data."""
        print(f"Loading data from {self.data_path}...")
        self.df = pd.read_csv(self.data_path)
        print(f"✓ Loaded {len(self.df)} data points\n")
        
    def run_symbolic_regression_full(self, 
                                     niterations=100,
                                     populations=30,
                                     population_size=50,
                                     maxsize=30,
                                     sample_frac=1.0):
        """
        Run symbolic regression on full dataset to discover governing equation.
        
        Parameters:
        -----------
        niterations : int
            Number of iterations (more = better but slower)
        populations : int
            Number of populations for evolution
        population_size : int
            Size of each population
        maxsize : int
            Maximum complexity of equations
        sample_frac : float
            Fraction of data to use (use <1.0 for faster iteration during testing)
        """
        if not PYSR_AVAILABLE:
            print("❌ PySR not installed!")
            return None
        
        print("="*70)
        print(" PHASE 3: SYMBOLIC REGRESSION - FULL DATASET")
        print("="*70)
        
        # Prepare features
        feature_cols = ['salinity', 'temperature', 'aou', 'stratification_index', 'c_atm']
        
        # Sample data if requested
        if sample_frac < 1.0:
            df_sample = self.df.sample(frac=sample_frac, random_state=42)
            print(f"\n📊 Using {len(df_sample)} samples ({sample_frac*100:.0f}% of data)")
        else:
            df_sample = self.df
            print(f"\n📊 Using full dataset: {len(df_sample)} samples")
        
        X = df_sample[feature_cols].values
        y = df_sample['co2_aq'].values
        
        # Train-test split
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )
        
        print(f"Train: {len(X_train)} samples")
        print(f"Test:  {len(X_test)} samples\n")
        
        # Feature scaling (helps PySR converge faster)
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Store for later
        self.X_train = X_train
        self.X_test = X_test
        self.y_train = y_train
        self.y_test = y_test
        self.feature_cols = feature_cols
        self.scaler = scaler
        
        print("🔬 Starting Symbolic Regression...")
        print(f"   Iterations: {niterations}")
        print(f"   Populations: {populations}")
        print(f"   Max complexity: {maxsize}")
        print(f"   Binary operators: [+, -, *, /]")
        print(f"   Unary operators: [exp, log, sqrt]\n")
        print("⏳ This may take 5-30 minutes depending on settings...\n")
        
        # Configure PySR
        model = PySRRegressor(
            niterations=niterations,
            populations=populations,
            population_size=population_size,
            binary_operators=["+", "-", "*", "/"],
            unary_operators=["exp", "log", "sqrt", "square"],
            maxsize=maxsize,
            constraints={
                '/': (-1, 5),
                'exp': 5,
                'log': 5
            },
            model_selection="best",
            verbosity=1,
            progress=True,
            random_state=42,
            temp_equation_file=True,
            parsimony=0.001
        )
        
        # Fit model with variable names
        print("="*70)
        model.fit(X_train_scaled, y_train, variable_names=feature_cols)
        print("="*70)
        
        # Get predictions
        y_pred_train = model.predict(X_train_scaled)
        y_pred_test = model.predict(X_test_scaled)
        
        # Metrics
        train_r2 = r2_score(y_train, y_pred_train)
        test_r2 = r2_score(y_test, y_pred_test)
        test_rmse = np.sqrt(mean_squared_error(y_test, y_pred_test))
        
        print("\n" + "="*70)
        print(" SYMBOLIC REGRESSION RESULTS")
        print("="*70)
        print(f"\n✓ Train R²:  {train_r2:.4f}")
        print(f"✓ Test R²:   {test_r2:.4f}")
        print(f"✓ Test RMSE: {test_rmse:.4f}")
        
        # Store results
        self.pysr_models['full'] = model
        self.results['full'] = {
            'model': model,
            'train_r2': train_r2,
            'test_r2': test_r2,
            'test_rmse': test_rmse,
            'predictions': y_pred_test
        }
        
        # Display equation hall of fame
        print("\n" + "="*70)
        print(" EQUATION HALL OF FAME (Top 10)")
        print("="*70 + "\n")
        
        try:
            equations = model.equations_
            print(equations[['complexity', 'loss', 'score', 'equation']].head(10).to_string(index=False))
        except Exception as e:
            print(f"(Unable to display equation table: {e})")
        
        # Get best equation
        print("\n" + "="*70)
        print(" 🏆 BEST EQUATION (Discovery)")
        print("="*70)
        
        try:
            best_eq = model.sympy()
            print(f"\n{best_eq}\n")
        except:
            try:
                best_eq = model.get_best()
                print(f"\n{best_eq}\n")
            except:
                print("\n(Equation display not available)\n")
        
        # Feature mapping
        print("\nFeature Mapping:")
        for i, name in enumerate(feature_cols):
            print(f"  x{i} = {name}")
        
        print("\n" + "="*70 + "\n")
        
        return model
    
    def run_symbolic_regression_simple(self, 
                                       niterations=80,
                                       maxsize=15):
        """
        Run simplified symbolic regression focusing on Salinity + key features.
        This produces simpler, more interpretable equations.
        """
        if not PYSR_AVAILABLE:
            print("❌ PySR not installed!")
            return None
        
        print("="*70)
        print(" PHASE 3B: SIMPLIFIED SYMBOLIC REGRESSION")
        print(" Focus: Salinity + Temperature + AOU")
        print("="*70)
        
        # Use only key features for simpler equation
        feature_cols_simple = ['salinity', 'temperature', 'aou']
        
        X = self.df[feature_cols_simple].values
        y = self.df['co2_aq'].values
        
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )
        
        print(f"\n📊 Using 3 core features:")
        for feat in feature_cols_simple:
            print(f"   - {feat}")
        print()
        
        # Scale
        scaler_simple = StandardScaler()
        X_train_scaled = scaler_simple.fit_transform(X_train)
        X_test_scaled = scaler_simple.transform(X_test)
        
        print("🔬 Running simplified search...")
        print(f"   Max complexity: {maxsize}")
        print(f"   Iterations: {niterations}\n")
        
        model_simple = PySRRegressor(
            niterations=niterations,
            populations=20,
            population_size=40,
            binary_operators=["+", "-", "*", "/"],
            unary_operators=["exp", "log", "sqrt"],
            maxsize=maxsize,
            model_selection="best",
            verbosity=1,
            progress=True,
            random_state=42,
            parsimony=0.01
        )
        
        # Fit with variable names
        model_simple.fit(X_train_scaled, y_train, variable_names=feature_cols_simple)
        
        # Evaluate
        y_pred_test = model_simple.predict(X_test_scaled)
        test_r2 = r2_score(y_test, y_pred_test)
        test_rmse = np.sqrt(mean_squared_error(y_test, y_pred_test))
        
        print("\n✓ Simple Model:")
        print(f"  Test R²:   {test_r2:.4f}")
        print(f"  Test RMSE: {test_rmse:.4f}")
        
        try:
            best_eq_simple = model_simple.sympy()
            print(f"\n🏆 Simplified Equation:\n\n{best_eq_simple}\n")
        except:
            try:
                best_eq_simple = model_simple.get_best()
                print(f"\n🏆 Simplified Equation:\n\n{best_eq_simple}\n")
            except:
                print("\n(Equation display not available)\n")
        
        print("Feature Mapping:")
        for i, name in enumerate(feature_cols_simple):
            print(f"  x{i} = {name}")
        
        self.pysr_models['simple'] = model_simple
        self.results['simple'] = {
            'model': model_simple,
            'test_r2': test_r2,
            'test_rmse': test_rmse,
            'features': feature_cols_simple,
            'scaler': scaler_simple,
            'predictions': y_pred_test
        }
        
        return model_simple
    
    def visualize_equation_performance(self):
        """
        Visualize how discovered equations perform.
        """
        print("\n" + "="*70)
        print(" Generating Visualizations")
        print("="*70 + "\n")
        
        model_names = []
        if 'full' in self.results:
            model_names.append('full')
        if 'simple' in self.results:
            model_names.append('simple')
        
        if len(model_names) == 0:
            print("⚠ No models to visualize")
            return
        
        fig, axes = plt.subplots(1, len(model_names), figsize=(7*len(model_names), 6))
        
        if len(model_names) == 1:
            axes = [axes]
        
        for idx, model_key in enumerate(model_names):
            result = self.results[model_key]
            y_pred = result['predictions']
            r2 = result['test_r2']
            
            # Get correct y_test
            if model_key == 'full':
                y_test = self.y_test
            else:
                y_test = self.df['co2_aq'].iloc[self.results['simple']['predictions'].index] if hasattr(self.results['simple']['predictions'], 'index') else self.y_test
            
            ax = axes[idx]
            
            # Scatter plot
            ax.scatter(y_test, y_pred, alpha=0.3, s=10, color='#3498db')
            
            # Perfect prediction line
            min_val = min(y_test.min(), y_pred.min())
            max_val = max(y_test.max(), y_pred.max())
            ax.plot([min_val, max_val], [min_val, max_val], 
                   'r--', linewidth=2, label='Perfect Prediction')
            
            label = 'Full Model (5 features)' if model_key == 'full' else 'Simple Model (3 features)'
            ax.set_xlabel('Actual CO₂(aq) (μmol/kg)', fontweight='bold')
            ax.set_ylabel('Predicted CO₂(aq) (μmol/kg)', fontweight='bold')
            ax.set_title(f'{label}\nR² = {r2:.4f}', fontweight='bold')
            ax.legend(loc='upper left')
            ax.grid(alpha=0.3)
        
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, '01_equation_performance.png')
        plt.savefig(save_path, bbox_inches='tight')
        print(f"✓ Saved: {save_path}")
        plt.close()
        
    def compare_with_ml_models(self):
        """
        Compare symbolic regression with ML models from Phase 2.
        """
        print("\nGenerating ML vs Symbolic comparison...")
        
        # Load Phase 2 results
        phase2_path = 'results/phase2/salinity_co2_ml_analysis/overall_results.csv'
        
        comparison_data = []
        
        if os.path.exists(phase2_path):
            df_phase2 = pd.read_csv(phase2_path)
            for _, row in df_phase2.iterrows():
                comparison_data.append({
                    'Model': row['model'],
                    'R²': row['r2'],
                    'RMSE': row['rmse'],
                    'Type': 'Black-box ML'
                })
        
        # Add symbolic regression results
        if 'full' in self.results:
            comparison_data.append({
                'Model': 'PySR (Full)',
                'R²': self.results['full']['test_r2'],
                'RMSE': self.results['full']['test_rmse'],
                'Type': 'Interpretable'
            })
        
        if 'simple' in self.results:
            comparison_data.append({
                'Model': 'PySR (Simple)',
                'R²': self.results['simple']['test_r2'],
                'RMSE': self.results['simple']['test_rmse'],
                'Type': 'Interpretable'
            })
        
        df_comp = pd.DataFrame(comparison_data)
        
        # Visualization
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # R² comparison
        ax1 = axes[0]
        colors = ['gray' if t == 'Black-box ML' else '#e74c3c' for t in df_comp['Type']]
        bars = ax1.bar(range(len(df_comp)), df_comp['R²'], color=colors, alpha=0.7)
        ax1.set_xticks(range(len(df_comp)))
        ax1.set_xticklabels(df_comp['Model'], rotation=45, ha='right')
        ax1.set_ylabel('R² Score', fontweight='bold')
        ax1.set_title('Model Comparison: R² Performance', fontweight='bold')
        ax1.axhline(y=0.99, color='green', linestyle='--', alpha=0.5, label='R²=0.99 threshold')
        ax1.legend()
        ax1.grid(axis='y', alpha=0.3)
        
        # Add value labels
        for i, bar in enumerate(bars):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.4f}',
                    ha='center', va='bottom', fontsize=8, fontweight='bold')
        
        # RMSE comparison
        ax2 = axes[1]
        bars = ax2.bar(range(len(df_comp)), df_comp['RMSE'], color=colors, alpha=0.7)
        ax2.set_xticks(range(len(df_comp)))
        ax2.set_xticklabels(df_comp['Model'], rotation=45, ha='right')
        ax2.set_ylabel('RMSE (μmol/kg)', fontweight='bold')
        ax2.set_title('Model Comparison: RMSE', fontweight='bold')
        ax2.grid(axis='y', alpha=0.3)
        
        # Add value labels
        for i, bar in enumerate(bars):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.2f}',
                    ha='center', va='bottom', fontsize=8, fontweight='bold')
        
        # Legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='gray', alpha=0.7, label='Black-box ML'),
            Patch(facecolor='#e74c3c', alpha=0.7, label='Interpretable (PySR)')
        ]
        ax2.legend(handles=legend_elements, loc='upper right')
        
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, '02_ml_vs_symbolic_comparison.png')
        plt.savefig(save_path, bbox_inches='tight')
        print(f"✓ Saved: {save_path}")
        plt.close()
        
        # Save comparison table
        csv_path = os.path.join(self.output_dir, 'model_comparison.csv')
        df_comp.to_csv(csv_path, index=False)
        print(f"✓ Saved: {csv_path}")
        
    def save_equations(self):
        """
        Save discovered equations in human-readable format.
        """
        print("\nSaving discovered equations...")
        
        eq_path = os.path.join(self.output_dir, 'DISCOVERED_EQUATIONS.txt')
        
        with open(eq_path, 'w') as f:
            f.write("="*70 + "\n")
            f.write(" DISCOVERED GOVERNING EQUATIONS\n")
            f.write(" Symbolic Regression with PySR\n")
            f.write("="*70 + "\n\n")
            
            if 'full' in self.results:
                f.write("## FULL MODEL (5 Features)\n")
                f.write("-"*70 + "\n\n")
                f.write(f"Features: salinity, temperature, aou, stratification_index, c_atm\n\n")
                f.write(f"Test R²: {self.results['full']['test_r2']:.4f}\n")
                f.write(f"Test RMSE: {self.results['full']['test_rmse']:.4f}\n\n")
                f.write("Equation:\n")
                try:
                    eq_str = str(self.pysr_models['full'].sympy())
                    f.write(f"{eq_str}\n\n")
                except:
                    try:
                        eq_str = str(self.pysr_models['full'].get_best())
                        f.write(f"{eq_str}\n\n")
                    except:
                        f.write("(Equation not available)\n\n")
                
                f.write("Feature Mapping:\n")
                for i, name in enumerate(self.feature_cols):
                    f.write(f"  x{i} = {name}\n")
                f.write("\n")
            
            if 'simple' in self.results:
                f.write("\n" + "="*70 + "\n\n")
                f.write("## SIMPLIFIED MODEL (3 Features)\n")
                f.write("-"*70 + "\n\n")
                f.write(f"Features: salinity, temperature, aou\n\n")
                f.write(f"Test R²: {self.results['simple']['test_r2']:.4f}\n")
                f.write(f"Test RMSE: {self.results['simple']['test_rmse']:.4f}\n\n")
                f.write("Equation:\n")
                try:
                    eq_str = str(self.pysr_models['simple'].sympy())
                    f.write(f"{eq_str}\n\n")
                except:
                    try:
                        eq_str = str(self.pysr_models['simple'].get_best())
                        f.write(f"{eq_str}\n\n")
                    except:
                        f.write("(Equation not available)\n\n")
                
                f.write("Feature Mapping:\n")
                for i, name in enumerate(self.results['simple']['features']):
                    f.write(f"  x{i} = {name}\n")
                f.write("\n")
            
            f.write("\n" + "="*70 + "\n")
            f.write(" SCIENTIFIC INTERPRETATION\n")
            f.write("="*70 + "\n\n")
            f.write("These equations represent the first data-driven, interpretable\n")
            f.write("mathematical model of the Salinity-CO₂(aq) relationship in the\n")
            f.write("Southern Ocean.\n\n")
            f.write("Key advantages over black-box ML:\n")
            f.write("1. Human-readable and verifiable\n")
            f.write("2. Can be integrated into numerical ocean models\n")
            f.write("3. Reveals physical mechanisms (e.g., exp/log terms)\n")
            f.write("4. Enables extrapolation beyond training data\n")
            f.write("5. Suitable for policy and climate impact assessments\n")
        
        print(f"✓ Saved: {eq_path}\n")
        
    def run_complete_phase3(self, 
                           quick_test=False,
                           run_simplified=True):
        """
        Execute complete Phase 3 workflow.
        
        Parameters:
        -----------
        quick_test : bool
            If True, use reduced settings for faster testing
        run_simplified : bool
            Whether to also run simplified 3-feature model
        """
        if not PYSR_AVAILABLE:
            print("\n" + "="*70)
            print(" ❌ PySR NOT INSTALLED")
            print("="*70)
            print("\nTo install PySR:")
            print("  1. pip install pysr")
            print("  2. python -m pysr install")
            print("\nThis will take 2-3 minutes.\n")
            return
        
        print("\n" + "="*70)
        print(" PHASE 3: SYMBOLIC REGRESSION")
        print(" Discovering Governing Equations")
        print("="*70 + "\n")
        
        self.load_data()
        
        # Set parameters based on mode
        if quick_test:
            print("⚡ QUICK TEST MODE (reduced settings)\n")
            self.run_symbolic_regression_full(
                niterations=30,
                populations=10,
                population_size=30,
                maxsize=20,
                sample_frac=0.3
            )
            if run_simplified:
                self.run_symbolic_regression_simple(
                    niterations=20,
                    maxsize=12
                )
        else:
            print("🔬 FULL ANALYSIS MODE\n")
            self.run_symbolic_regression_full(
                niterations=100,
                populations=30,
                population_size=50,
                maxsize=30,
                sample_frac=1.0
            )
            if run_simplified:
                self.run_symbolic_regression_simple(
                    niterations=80,
                    maxsize=15
                )
        
        # Generate outputs
        self.visualize_equation_performance()
        self.compare_with_ml_models()
        self.save_equations()
        
        print("\n" + "="*70)
        print(" ✓ PHASE 3 COMPLETE!")
        print("="*70)
        print(f"\nResults saved to: {self.output_dir}/")
        print("\nGenerated files:")
        print("  01_equation_performance.png       - Discovered equation performance")
        print("  02_ml_vs_symbolic_comparison.png  - ML vs Symbolic comparison")
        print("  DISCOVERED_EQUATIONS.txt          - Human-readable equations")
        print("  model_comparison.csv              - Numerical comparison")
        print("\n🎯 You now have an interpretable governing equation!")
        print("   Ready for Phase 4: Physical interpretation\n")


def main():
    """
    Main execution.
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Phase 3: Symbolic Regression")
    parser.add_argument('--quick', action='store_true', 
                       help='Quick test mode (30 iterations, 30%% data)')
    parser.add_argument('--full-only', action='store_true',
                       help='Skip simplified 3-feature model')
    
    args = parser.parse_args()
    
    phase3 = SymbolicRegressionPhase3(
        data_path='data/processed/southern_ocean_training.csv',
        output_dir='results/phase3'
    )
    
    phase3.run_complete_phase3(
        quick_test=args.quick,
        run_simplified=not args.full_only
    )


if __name__ == "__main__":
    main()
