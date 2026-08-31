"""
Baseline Machine Learning Models for Phase 2
- Random Forest Regressor
- XGBoost Regressor
- Feature Importance Analysis with SHAP
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import joblib
import os

try:
    import xgboost as xgb
    XGBOOST_AVAILABLE = True
except ImportError:
    XGBOOST_AVAILABLE = False
    print("Warning: XGBoost not installed. Install with: pip install xgboost")

try:
    import shap
    SHAP_AVAILABLE = True
except ImportError:
    SHAP_AVAILABLE = False
    print("Warning: SHAP not installed. Install with: pip install shap")


class BaselineModelTrainer:
    """
    Trains and evaluates baseline regression models for CO2(aq) prediction.
    """
    
    def __init__(self, data_path, target_col='co2_aq', test_size=0.2, random_state=42):
        """
        Initialize the baseline trainer.
        
        Parameters:
        -----------
        data_path : str
            Path to the processed CSV file
        target_col : str
            Name of target variable column
        test_size : float
            Fraction of data for testing
        random_state : int
            Random seed for reproducibility
        """
        self.data_path = data_path
        self.target_col = target_col
        self.test_size = test_size
        self.random_state = random_state
        
        # Data containers
        self.df = None
        self.X_train = None
        self.X_test = None
        self.y_train = None
        self.y_test = None
        self.feature_names = None
        
        # Model containers
        self.models = {}
        self.results = {}
        
    def load_and_prepare_data(self, exclude_cols=None):
        """
        Load data and split into train/test sets.
        
        Parameters:
        -----------
        exclude_cols : list
            Columns to exclude from features (beyond target)
        """
        print(f"Loading data from {self.data_path}...")
        self.df = pd.read_csv(self.data_path)
        print(f"Dataset loaded: {len(self.df)} rows, {len(self.df.columns)} columns")
        
        # Define feature columns
        if exclude_cols is None:
            exclude_cols = ['year', 'latitude', 'longitude', 'depth', 'pressure']
        
        # Features: Everything except target and excluded columns
        feature_cols = [col for col in self.df.columns 
                       if col != self.target_col and col not in exclude_cols]
        
        X = self.df[feature_cols]
        y = self.df[self.target_col]
        
        self.feature_names = feature_cols
        
        # Train-test split
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(
            X, y, test_size=self.test_size, random_state=self.random_state
        )
        
        print(f"\nFeatures used ({len(feature_cols)}):")
        for i, col in enumerate(feature_cols, 1):
            print(f"  {i}. {col}")
        
        print(f"\nTrain set: {len(self.X_train)} samples")
        print(f"Test set:  {len(self.X_test)} samples")
        
    def train_random_forest(self, n_estimators=200, max_depth=None, 
                           min_samples_split=5, n_jobs=-1):
        """
        Train Random Forest model.
        
        Parameters:
        -----------
        n_estimators : int
            Number of trees
        max_depth : int or None
            Maximum tree depth
        min_samples_split : int
            Minimum samples required to split
        n_jobs : int
            Number of parallel jobs (-1 = all cores)
        """
        print("\n" + "="*60)
        print("Training Random Forest Regressor...")
        print("="*60)
        
        rf = RandomForestRegressor(
            n_estimators=n_estimators,
            max_depth=max_depth,
            min_samples_split=min_samples_split,
            random_state=self.random_state,
            n_jobs=n_jobs,
            verbose=0
        )
        
        rf.fit(self.X_train, self.y_train)
        
        # Predictions
        y_train_pred = rf.predict(self.X_train)
        y_test_pred = rf.predict(self.X_test)
        
        # Metrics
        train_r2 = r2_score(self.y_train, y_train_pred)
        test_r2 = r2_score(self.y_test, y_test_pred)
        test_rmse = np.sqrt(mean_squared_error(self.y_test, y_test_pred))
        test_mae = mean_absolute_error(self.y_test, y_test_pred)
        
        # Feature importance
        feature_importance = pd.DataFrame({
            'feature': self.feature_names,
            'importance': rf.feature_importances_
        }).sort_values('importance', ascending=False)
        
        # Store results
        self.models['RandomForest'] = rf
        self.results['RandomForest'] = {
            'train_r2': train_r2,
            'test_r2': test_r2,
            'test_rmse': test_rmse,
            'test_mae': test_mae,
            'feature_importance': feature_importance,
            'predictions': y_test_pred
        }
        
        print(f"\n✓ Random Forest Results:")
        print(f"  Train R²:  {train_r2:.4f}")
        print(f"  Test R²:   {test_r2:.4f}")
        print(f"  Test RMSE: {test_rmse:.4f}")
        print(f"  Test MAE:  {test_mae:.4f}")
        
        print(f"\nTop 5 Most Important Features:")
        for idx, row in feature_importance.head(5).iterrows():
            print(f"  {row['feature']:25s} {row['importance']:.4f}")
            
        return rf, test_r2
    
    def train_xgboost(self, n_estimators=200, max_depth=6, 
                     learning_rate=0.1, subsample=0.8):
        """
        Train XGBoost model.
        
        Parameters:
        -----------
        n_estimators : int
            Number of boosting rounds
        max_depth : int
            Maximum tree depth
        learning_rate : float
            Step size shrinkage
        subsample : float
            Subsample ratio of training data
        """
        if not XGBOOST_AVAILABLE:
            print("\n⚠ XGBoost not available. Skipping.")
            return None, None
        
        print("\n" + "="*60)
        print("Training XGBoost Regressor...")
        print("="*60)
        
        xgb_model = xgb.XGBRegressor(
            n_estimators=n_estimators,
            max_depth=max_depth,
            learning_rate=learning_rate,
            subsample=subsample,
            random_state=self.random_state,
            n_jobs=-1,
            verbosity=0
        )
        
        xgb_model.fit(self.X_train, self.y_train)
        
        # Predictions
        y_train_pred = xgb_model.predict(self.X_train)
        y_test_pred = xgb_model.predict(self.X_test)
        
        # Metrics
        train_r2 = r2_score(self.y_train, y_train_pred)
        test_r2 = r2_score(self.y_test, y_test_pred)
        test_rmse = np.sqrt(mean_squared_error(self.y_test, y_test_pred))
        test_mae = mean_absolute_error(self.y_test, y_test_pred)
        
        # Feature importance
        feature_importance = pd.DataFrame({
            'feature': self.feature_names,
            'importance': xgb_model.feature_importances_
        }).sort_values('importance', ascending=False)
        
        # Store results
        self.models['XGBoost'] = xgb_model
        self.results['XGBoost'] = {
            'train_r2': train_r2,
            'test_r2': test_r2,
            'test_rmse': test_rmse,
            'test_mae': test_mae,
            'feature_importance': feature_importance,
            'predictions': y_test_pred
        }
        
        print(f"\n✓ XGBoost Results:")
        print(f"  Train R²:  {train_r2:.4f}")
        print(f"  Test R²:   {test_r2:.4f}")
        print(f"  Test RMSE: {test_rmse:.4f}")
        print(f"  Test MAE:  {test_mae:.4f}")
        
        print(f"\nTop 5 Most Important Features:")
        for idx, row in feature_importance.head(5).iterrows():
            print(f"  {row['feature']:25s} {row['importance']:.4f}")
            
        return xgb_model, test_r2
    
    def compute_shap_values(self, model_name='RandomForest', sample_size=500):
        """
        Compute SHAP values for model interpretability.
        
        Parameters:
        -----------
        model_name : str
            Which model to analyze
        sample_size : int
            Number of samples for SHAP (computational efficiency)
        """
        if not SHAP_AVAILABLE:
            print("\n⚠ SHAP not available. Skipping.")
            return None
        
        if model_name not in self.models:
            print(f"\n⚠ Model '{model_name}' not trained yet.")
            return None
        
        print(f"\n" + "="*60)
        print(f"Computing SHAP values for {model_name}...")
        print("="*60)
        
        model = self.models[model_name]
        
        # Sample data for efficiency
        if len(self.X_test) > sample_size:
            sample_idx = np.random.choice(len(self.X_test), sample_size, replace=False)
            X_sample = self.X_test.iloc[sample_idx]
        else:
            X_sample = self.X_test
        
        # Create explainer
        if model_name == 'XGBoost' and XGBOOST_AVAILABLE:
            explainer = shap.TreeExplainer(model)
        else:
            explainer = shap.TreeExplainer(model)
        
        shap_values = explainer.shap_values(X_sample)
        
        # Store SHAP results
        self.results[model_name]['shap_values'] = shap_values
        self.results[model_name]['shap_data'] = X_sample
        
        # Summary statistics
        mean_abs_shap = np.abs(shap_values).mean(axis=0)
        shap_importance = pd.DataFrame({
            'feature': self.feature_names,
            'mean_abs_shap': mean_abs_shap
        }).sort_values('mean_abs_shap', ascending=False)
        
        print(f"\n✓ SHAP Analysis Complete")
        print(f"\nTop 5 Features by SHAP Impact:")
        for idx, row in shap_importance.head(5).iterrows():
            print(f"  {row['feature']:25s} {row['mean_abs_shap']:.4f}")
        
        return shap_values
    
    def save_models(self, save_dir='models/phase2'):
        """
        Save trained models to disk.
        
        Parameters:
        -----------
        save_dir : str
            Directory to save models
        """
        os.makedirs(save_dir, exist_ok=True)
        
        for name, model in self.models.items():
            filepath = os.path.join(save_dir, f"{name.lower()}_model.joblib")
            joblib.dump(model, filepath)
            print(f"✓ Saved {name} to {filepath}")
    
    def save_results(self, save_dir='results/phase2'):
        """
        Save results and metrics to CSV files.
        
        Parameters:
        -----------
        save_dir : str
            Directory to save results
        """
        os.makedirs(save_dir, exist_ok=True)
        
        # Summary metrics
        summary = []
        for name, res in self.results.items():
            summary.append({
                'model': name,
                'train_r2': res['train_r2'],
                'test_r2': res['test_r2'],
                'test_rmse': res['test_rmse'],
                'test_mae': res['test_mae']
            })
        
        summary_df = pd.DataFrame(summary)
        summary_path = os.path.join(save_dir, 'model_comparison.csv')
        summary_df.to_csv(summary_path, index=False)
        print(f"\n✓ Saved model comparison to {summary_path}")
        
        # Feature importance for each model
        for name, res in self.results.items():
            if 'feature_importance' in res:
                fi_path = os.path.join(save_dir, f'{name.lower()}_feature_importance.csv')
                res['feature_importance'].to_csv(fi_path, index=False)
                print(f"✓ Saved {name} feature importance to {fi_path}")


if __name__ == "__main__":
    # Example usage
    trainer = BaselineModelTrainer(
        data_path='data/processed/southern_ocean_training.csv',
        test_size=0.2,
        random_state=42
    )
    
    trainer.load_and_prepare_data()
    trainer.train_random_forest(n_estimators=200)
    trainer.train_xgboost(n_estimators=200)
    trainer.compute_shap_values('RandomForest', sample_size=500)
    
    trainer.save_models()
    trainer.save_results()
