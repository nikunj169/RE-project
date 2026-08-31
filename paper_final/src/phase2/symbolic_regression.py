"""
Symbolic Regression using PySR for discovering governing equations
"""

import numpy as np
import pandas as pd
import os

try:
    from pysr import PySRRegressor
    PYSR_AVAILABLE = True
except ImportError:
    PYSR_AVAILABLE = False
    print("Warning: PySR not installed. Install with: pip install pysr")


class SymbolicRegressionDiscovery:
    """
    Discover symbolic equations for CO2(aq) using PySR.
    """
    
    def __init__(self, data_path, target_col='co2_aq', test_size=0.2, random_state=42):
        """
        Initialize symbolic regression.
        
        Parameters:
        -----------
        data_path : str
            Path to processed CSV
        target_col : str
            Target variable name
        test_size : float
            Test set fraction
        random_state : int
            Random seed
        """
        self.data_path = data_path
        self.target_col = target_col
        self.test_size = test_size
        self.random_state = random_state
        
        self.df = None
        self.X_train = None
        self.X_test = None
        self.y_train = None
        self.y_test = None
        self.feature_names = None
        
        self.pysr_model = None
        self.best_equation = None
        
    def load_and_prepare_data(self, exclude_cols=None, sample_frac=1.0):
        """
        Load and prepare data for symbolic regression.
        
        Parameters:
        -----------
        exclude_cols : list
            Columns to exclude from features
        sample_frac : float
            Fraction of data to use (for computational efficiency)
        """
        print(f"Loading data from {self.data_path}...")
        self.df = pd.read_csv(self.data_path)
        
        # Optional sampling for faster iteration
        if sample_frac < 1.0:
            self.df = self.df.sample(frac=sample_frac, random_state=self.random_state)
            print(f"Sampled {len(self.df)} rows ({sample_frac*100:.1f}% of data)")
        
        # Define features
        if exclude_cols is None:
            exclude_cols = ['year', 'latitude', 'longitude', 'depth', 'pressure']
        
        feature_cols = [col for col in self.df.columns 
                       if col != self.target_col and col not in exclude_cols]
        
        X = self.df[feature_cols]
        y = self.df[self.target_col]
        
        self.feature_names = feature_cols
        
        # Train-test split
        from sklearn.model_selection import train_test_split
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(
            X, y, test_size=self.test_size, random_state=self.random_state
        )
        
        print(f"\nFeatures for symbolic regression ({len(feature_cols)}):")
        for i, col in enumerate(feature_cols, 1):
            print(f"  {i}. {col}")
        
        print(f"\nTrain set: {len(self.X_train)} samples")
        print(f"Test set:  {len(self.X_test)} samples")
        
    def run_symbolic_regression(self, 
                                populations=15,
                                population_size=50,
                                niterations=50,
                                maxsize=20,
                                binary_operators=None,
                                unary_operators=None):
        """
        Run PySR to discover symbolic equation.
        
        Parameters:
        -----------
        populations : int
            Number of populations
        population_size : int
            Population size per generation
        niterations : int
            Number of iterations
        maxsize : int
            Maximum complexity of equations
        binary_operators : list
            Allowed binary operators
        unary_operators : list
            Allowed unary operators
        """
        if not PYSR_AVAILABLE:
            print("\n⚠ PySR not installed. Please install:")
            print("   pip install pysr")
            print("   python -m pysr install")
            return None
        
        if binary_operators is None:
            # Exclude sin/cos as they rarely appear in chemical thermodynamics
            binary_operators = ["+", "-", "*", "/"]
        
        if unary_operators is None:
            # exp, log are common in oceanographic equations
            unary_operators = ["exp", "log", "sqrt"]
        
        print("\n" + "="*60)
        print("Running PySR Symbolic Regression...")
        print("="*60)
        print(f"Binary operators: {binary_operators}")
        print(f"Unary operators:  {unary_operators}")
        print(f"Max complexity:   {maxsize}")
        print(f"Iterations:       {niterations}")
        print("\nThis may take several minutes...\n")
        
        self.pysr_model = PySRRegressor(
            populations=populations,
            population_size=population_size,
            niterations=niterations,
            binary_operators=binary_operators,
            unary_operators=unary_operators,
            maxsize=maxsize,
            model_selection="best",
            loss="loss(x, y) = (x - y)^2",
            random_state=self.random_state,
            temp_equation_file=True,
            verbosity=1,
            progress=True
        )
        
        # Fit the model
        self.pysr_model.fit(self.X_train.values, self.y_train.values)
        
        # Get best equation
        print("\n" + "="*60)
        print("PySR Search Complete!")
        print("="*60)
        
        # Display equation hall of fame
        print("\nEquation Hall of Fame (sorted by complexity vs accuracy):")
        print(self.pysr_model.equations_)
        
        # Get best equation
        self.best_equation = self.pysr_model.get_best()
        
        print("\n" + "="*60)
        print("BEST EQUATION:")
        print("="*60)
        print(f"\n{self.best_equation}\n")
        
        # Evaluate on test set
        y_pred_train = self.pysr_model.predict(self.X_train.values)
        y_pred_test = self.pysr_model.predict(self.X_test.values)
        
        from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
        
        train_r2 = r2_score(self.y_train, y_pred_train)
        test_r2 = r2_score(self.y_test, y_pred_test)
        test_rmse = np.sqrt(mean_squared_error(self.y_test, y_pred_test))
        test_mae = mean_absolute_error(self.y_test, y_pred_test)
        
        print(f"✓ Symbolic Regression Performance:")
        print(f"  Train R²:  {train_r2:.4f}")
        print(f"  Test R²:   {test_r2:.4f}")
        print(f"  Test RMSE: {test_rmse:.4f}")
        print(f"  Test MAE:  {test_mae:.4f}")
        
        self.results = {
            'train_r2': train_r2,
            'test_r2': test_r2,
            'test_rmse': test_rmse,
            'test_mae': test_mae,
            'best_equation': str(self.best_equation),
            'predictions_test': y_pred_test
        }
        
        return self.pysr_model
    
    def save_model(self, save_dir='models/phase2'):
        """
        Save PySR model.
        
        Parameters:
        -----------
        save_dir : str
            Directory to save model
        """
        if self.pysr_model is None:
            print("⚠ No model to save. Run symbolic regression first.")
            return
        
        os.makedirs(save_dir, exist_ok=True)
        filepath = os.path.join(save_dir, 'pysr_model.pkl')
        
        import pickle
        with open(filepath, 'wb') as f:
            pickle.dump(self.pysr_model, f)
        
        print(f"✓ Saved PySR model to {filepath}")
        
        # Also save equation as text
        eq_path = os.path.join(save_dir, 'discovered_equation.txt')
        with open(eq_path, 'w') as f:
            f.write("DISCOVERED EQUATION FOR CO2(aq)\n")
            f.write("="*60 + "\n\n")
            f.write(str(self.best_equation) + "\n\n")
            f.write("FEATURE MAPPING:\n")
            for i, name in enumerate(self.feature_names):
                f.write(f"  x{i} = {name}\n")
        
        print(f"✓ Saved equation to {eq_path}")
    
    def save_results(self, save_dir='results/phase2'):
        """
        Save symbolic regression results.
        
        Parameters:
        -----------
        save_dir : str
            Directory to save results
        """
        if self.results is None:
            print("⚠ No results to save. Run symbolic regression first.")
            return
        
        os.makedirs(save_dir, exist_ok=True)
        
        # Save metrics
        metrics_df = pd.DataFrame([{
            'model': 'PySR_SymbolicRegression',
            'train_r2': self.results['train_r2'],
            'test_r2': self.results['test_r2'],
            'test_rmse': self.results['test_rmse'],
            'test_mae': self.results['test_mae'],
            'equation': self.results['best_equation']
        }])
        
        metrics_path = os.path.join(save_dir, 'pysr_results.csv')
        metrics_df.to_csv(metrics_path, index=False)
        print(f"✓ Saved PySR metrics to {metrics_path}")


if __name__ == "__main__":
    # Example usage
    sr = SymbolicRegressionDiscovery(
        data_path='data/processed/southern_ocean_training.csv',
        test_size=0.2,
        random_state=42
    )
    
    # For faster iteration during development, use sample_frac < 1.0
    sr.load_and_prepare_data(sample_frac=0.3)
    
    sr.run_symbolic_regression(
        populations=10,
        population_size=30,
        niterations=30,
        maxsize=15
    )
    
    sr.save_model()
    sr.save_results()
