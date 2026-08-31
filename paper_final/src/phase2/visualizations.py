"""
Visualization suite for Phase 2 results
Generates publication-quality plots for conference paper
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
import os
from sklearn.metrics import r2_score
import warnings
warnings.filterwarnings('ignore')

# Set publication-quality style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9


class Phase2Visualizer:
    """
    Generate all visualizations for Phase 2 analysis.
    """
    
    def __init__(self, 
                 data_path='data/processed/southern_ocean_training.csv',
                 results_dir='results/phase2',
                 models_dir='models/phase2',
                 output_dir='results/phase2/figures'):
        """
        Initialize visualizer.
        
        Parameters:
        -----------
        data_path : str
            Path to processed training data
        results_dir : str
            Directory containing result CSVs
        models_dir : str
            Directory containing trained models
        output_dir : str
            Directory to save figures
        """
        self.data_path = data_path
        self.results_dir = results_dir
        self.models_dir = models_dir
        self.output_dir = output_dir
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Data containers
        self.df = None
        self.X_test = None
        self.y_test = None
        self.models = {}
        self.predictions = {}
        
    def load_data_and_models(self):
        """
        Load dataset and trained models.
        """
        print("Loading data and models...")
        
        # Load data
        self.df = pd.read_csv(self.data_path)
        
        # Prepare test set (same split as training)
        from sklearn.model_selection import train_test_split
        
        exclude_cols = ['year', 'latitude', 'longitude', 'depth', 'pressure', 'co2_aq']
        feature_cols = [col for col in self.df.columns if col not in exclude_cols]
        
        X = self.df[feature_cols]
        y = self.df['co2_aq']
        
        _, self.X_test, _, self.y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )
        
        # Load models
        rf_path = os.path.join(self.models_dir, 'randomforest_model.joblib')
        xgb_path = os.path.join(self.models_dir, 'xgboost_model.joblib')
        
        if os.path.exists(rf_path):
            self.models['Random Forest'] = joblib.load(rf_path)
            self.predictions['Random Forest'] = self.models['Random Forest'].predict(self.X_test)
            print(f"✓ Loaded Random Forest model")
        
        if os.path.exists(xgb_path):
            self.models['XGBoost'] = joblib.load(xgb_path)
            self.predictions['XGBoost'] = self.models['XGBoost'].predict(self.X_test)
            print(f"✓ Loaded XGBoost model")
        
        print(f"✓ Test set: {len(self.y_test)} samples\n")
    
    def plot_model_comparison(self):
        """
        Plot 1: Model Performance Comparison (Bar Chart)
        """
        print("Generating Plot 1: Model Performance Comparison...")
        
        comparison_path = os.path.join(self.results_dir, 'model_comparison.csv')
        df_comp = pd.read_csv(comparison_path)
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # R² Comparison
        ax1 = axes[0]
        x_pos = np.arange(len(df_comp))
        bars1 = ax1.bar(x_pos, df_comp['test_r2'], color=['#3498db', '#e74c3c'], alpha=0.8)
        ax1.set_ylabel('R² Score', fontweight='bold')
        ax1.set_title('Model Performance: R² on Test Set', fontweight='bold')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels(df_comp['model'], rotation=0)
        ax1.set_ylim([0.98, 1.0])
        ax1.axhline(y=0.99, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        ax1.grid(axis='y', alpha=0.3)
        
        # Add value labels on bars
        for bar in bars1:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.4f}',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        # RMSE Comparison
        ax2 = axes[1]
        bars2 = ax2.bar(x_pos, df_comp['test_rmse'], color=['#3498db', '#e74c3c'], alpha=0.8)
        ax2.set_ylabel('RMSE (μmol/kg)', fontweight='bold')
        ax2.set_title('Model Performance: RMSE on Test Set', fontweight='bold')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(df_comp['model'], rotation=0)
        ax2.grid(axis='y', alpha=0.3)
        
        # Add value labels on bars
        for bar in bars2:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.2f}',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, '01_model_comparison.png')
        plt.savefig(save_path, bbox_inches='tight')
        print(f"  ✓ Saved to {save_path}")
        plt.close()
    
    def plot_predicted_vs_actual(self):
        """
        Plot 2: Predicted vs Actual CO2(aq) for both models
        """
        print("Generating Plot 2: Predicted vs Actual...")
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        for idx, (model_name, y_pred) in enumerate(self.predictions.items()):
            ax = axes[idx]
            
            # Scatter plot
            scatter = ax.scatter(self.y_test, y_pred, 
                               alpha=0.3, s=10, c=self.y_test, cmap='viridis')
            
            # Perfect prediction line
            min_val = min(self.y_test.min(), y_pred.min())
            max_val = max(self.y_test.max(), y_pred.max())
            ax.plot([min_val, max_val], [min_val, max_val], 
                   'r--', linewidth=2, label='Perfect Prediction')
            
            # Calculate R²
            r2 = r2_score(self.y_test, y_pred)
            
            ax.set_xlabel('Actual CO₂(aq) (μmol/kg)', fontweight='bold')
            ax.set_ylabel('Predicted CO₂(aq) (μmol/kg)', fontweight='bold')
            ax.set_title(f'{model_name}: R² = {r2:.4f}', fontweight='bold')
            ax.legend(loc='upper left')
            ax.grid(alpha=0.3)
            
            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('Actual CO₂(aq)', rotation=270, labelpad=20)
        
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, '02_predicted_vs_actual.png')
        plt.savefig(save_path, bbox_inches='tight')
        print(f"  ✓ Saved to {save_path}")
        plt.close()
    
    def plot_residual_analysis(self):
        """
        Plot 3: Residual Analysis for both models
        """
        print("Generating Plot 3: Residual Analysis...")
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        for idx, (model_name, y_pred) in enumerate(self.predictions.items()):
            residuals = self.y_test.values - y_pred
            
            # Residuals vs Predicted (top row)
            ax1 = axes[0, idx]
            ax1.scatter(y_pred, residuals, alpha=0.3, s=10)
            ax1.axhline(y=0, color='r', linestyle='--', linewidth=2)
            ax1.set_xlabel('Predicted CO₂(aq) (μmol/kg)', fontweight='bold')
            ax1.set_ylabel('Residuals (μmol/kg)', fontweight='bold')
            ax1.set_title(f'{model_name}: Residual Plot', fontweight='bold')
            ax1.grid(alpha=0.3)
            
            # Residual histogram (bottom row)
            ax2 = axes[1, idx]
            ax2.hist(residuals, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
            ax2.axvline(x=0, color='r', linestyle='--', linewidth=2)
            ax2.set_xlabel('Residuals (μmol/kg)', fontweight='bold')
            ax2.set_ylabel('Frequency', fontweight='bold')
            ax2.set_title(f'{model_name}: Residual Distribution', fontweight='bold')
            ax2.grid(alpha=0.3)
            
            # Add statistics
            mean_res = np.mean(residuals)
            std_res = np.std(residuals)
            ax2.text(0.98, 0.95, f'Mean: {mean_res:.3f}\nStd: {std_res:.3f}',
                    transform=ax2.transAxes, ha='right', va='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, '03_residual_analysis.png')
        plt.savefig(save_path, bbox_inches='tight')
        print(f"  ✓ Saved to {save_path}")
        plt.close()
    
    def plot_feature_importance(self):
        """
        Plot 4: Feature Importance Comparison
        """
        print("Generating Plot 4: Feature Importance...")
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Random Forest
        rf_fi_path = os.path.join(self.results_dir, 'randomforest_feature_importance.csv')
        if os.path.exists(rf_fi_path):
            df_rf = pd.read_csv(rf_fi_path).head(8)
            
            ax1 = axes[0]
            bars1 = ax1.barh(df_rf['feature'], df_rf['importance'], color='#3498db', alpha=0.8)
            ax1.set_xlabel('Importance Score', fontweight='bold')
            ax1.set_title('Random Forest: Feature Importance', fontweight='bold')
            ax1.invert_yaxis()
            ax1.grid(axis='x', alpha=0.3)
            
            # Add value labels
            for i, bar in enumerate(bars1):
                width = bar.get_width()
                ax1.text(width, bar.get_y() + bar.get_height()/2.,
                        f'{width:.4f}',
                        ha='left', va='center', fontsize=9)
        
        # XGBoost
        xgb_fi_path = os.path.join(self.results_dir, 'xgboost_feature_importance.csv')
        if os.path.exists(xgb_fi_path):
            df_xgb = pd.read_csv(xgb_fi_path).head(8)
            
            ax2 = axes[1]
            bars2 = ax2.barh(df_xgb['feature'], df_xgb['importance'], color='#e74c3c', alpha=0.8)
            ax2.set_xlabel('Importance Score', fontweight='bold')
            ax2.set_title('XGBoost: Feature Importance', fontweight='bold')
            ax2.invert_yaxis()
            ax2.grid(axis='x', alpha=0.3)
            
            # Add value labels
            for i, bar in enumerate(bars2):
                width = bar.get_width()
                ax2.text(width, bar.get_y() + bar.get_height()/2.,
                        f'{width:.4f}',
                        ha='left', va='center', fontsize=9)
        
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, '04_feature_importance.png')
        plt.savefig(save_path, bbox_inches='tight')
        print(f"  ✓ Saved to {save_path}")
        plt.close()
    
    def plot_error_by_depth(self):
        """
        Plot 5: Model Error Analysis by Depth
        """
        print("Generating Plot 5: Error by Depth...")
        
        # Get depth information for test set
        _, X_test_full, _, y_test_full = train_test_split(
            self.df, self.df['co2_aq'], test_size=0.2, random_state=42
        )
        
        depths = X_test_full['depth'].values
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        for idx, (model_name, y_pred) in enumerate(self.predictions.items()):
            ax = axes[idx]
            
            # Calculate absolute errors
            abs_errors = np.abs(self.y_test.values - y_pred)
            
            # Bin by depth
            depth_bins = [0, 100, 500, 1000, 2000, 3000, 6000]
            depth_labels = ['0-100m', '100-500m', '500-1000m', 
                           '1000-2000m', '2000-3000m', '3000-6000m']
            
            depth_binned = pd.cut(depths, bins=depth_bins, labels=depth_labels)
            
            df_errors = pd.DataFrame({
                'depth_bin': depth_binned,
                'abs_error': abs_errors
            })
            
            # Box plot
            df_errors.boxplot(column='abs_error', by='depth_bin', ax=ax, patch_artist=True)
            ax.set_xlabel('Depth Range', fontweight='bold')
            ax.set_ylabel('Absolute Error (μmol/kg)', fontweight='bold')
            ax.set_title(f'{model_name}: Prediction Error by Depth', fontweight='bold')
            plt.sca(ax)
            plt.xticks(rotation=45, ha='right')
            ax.get_figure().suptitle('')  # Remove automatic title
        
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, '05_error_by_depth.png')
        plt.savefig(save_path, bbox_inches='tight')
        print(f"  ✓ Saved to {save_path}")
        plt.close()
    
    def plot_error_by_sector(self):
        """
        Plot 6: Model Error Analysis by Ocean Sector
        """
        print("Generating Plot 6: Error by Sector...")
        
        from sklearn.model_selection import train_test_split
        
        # Get sector information for test set
        _, X_test_full, _, y_test_full = train_test_split(
            self.df, self.df['co2_aq'], test_size=0.2, random_state=42
        )
        
        # Determine sector
        sectors = []
        for _, row in X_test_full.iterrows():
            if row['sector_atlantic'] == 1:
                sectors.append('Atlantic')
            elif row['sector_indian'] == 1:
                sectors.append('Indian')
            else:
                sectors.append('Pacific')
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        for idx, (model_name, y_pred) in enumerate(self.predictions.items()):
            ax = axes[idx]
            
            # Calculate errors
            abs_errors = np.abs(self.y_test.values - y_pred)
            
            df_errors = pd.DataFrame({
                'sector': sectors,
                'abs_error': abs_errors
            })
            
            # Box plot
            df_errors.boxplot(column='abs_error', by='sector', ax=ax, patch_artist=True)
            ax.set_xlabel('Ocean Sector', fontweight='bold')
            ax.set_ylabel('Absolute Error (μmol/kg)', fontweight='bold')
            ax.set_title(f'{model_name}: Prediction Error by Sector', fontweight='bold')
            ax.get_figure().suptitle('')  # Remove automatic title
            
            # Add sample sizes
            sector_counts = df_errors['sector'].value_counts()
            labels = [f'{s}\n(n={sector_counts[s]})' for s in ['Atlantic', 'Indian', 'Pacific']]
            ax.set_xticklabels(labels)
        
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, '06_error_by_sector.png')
        plt.savefig(save_path, bbox_inches='tight')
        print(f"  ✓ Saved to {save_path}")
        plt.close()
    
    def plot_feature_correlation_heatmap(self):
        """
        Plot 7: Feature Correlation Heatmap
        """
        print("Generating Plot 7: Feature Correlation Heatmap...")
        
        # Select numeric features
        exclude_cols = ['year', 'latitude', 'longitude', 'depth', 'pressure']
        feature_cols = [col for col in self.df.columns if col not in exclude_cols]
        
        df_features = self.df[feature_cols]
        
        # Calculate correlation matrix
        corr_matrix = df_features.corr()
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create heatmap
        sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm', 
                   center=0, square=True, linewidths=1, cbar_kws={"shrink": 0.8},
                   ax=ax)
        
        ax.set_title('Feature Correlation Matrix', fontweight='bold', fontsize=14, pad=20)
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, '07_correlation_heatmap.png')
        plt.savefig(save_path, bbox_inches='tight')
        print(f"  ✓ Saved to {save_path}")
        plt.close()
    
    def plot_learning_curves(self):
        """
        Plot 8: Prediction Distribution Comparison
        """
        print("Generating Plot 8: Prediction Distribution...")
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Plot actual distribution
        ax.hist(self.y_test, bins=50, alpha=0.5, label='Actual', 
               color='gray', edgecolor='black')
        
        # Plot predicted distributions
        colors = ['#3498db', '#e74c3c']
        for idx, (model_name, y_pred) in enumerate(self.predictions.items()):
            ax.hist(y_pred, bins=50, alpha=0.5, label=f'{model_name} Predicted',
                   color=colors[idx], edgecolor='black')
        
        ax.set_xlabel('CO₂(aq) (μmol/kg)', fontweight='bold')
        ax.set_ylabel('Frequency', fontweight='bold')
        ax.set_title('Distribution of Actual vs Predicted CO₂(aq)', fontweight='bold')
        ax.legend(loc='upper right')
        ax.grid(alpha=0.3)
        
        plt.tight_layout()
        save_path = os.path.join(self.output_dir, '08_prediction_distribution.png')
        plt.savefig(save_path, bbox_inches='tight')
        print(f"  ✓ Saved to {save_path}")
        plt.close()
    
    def generate_all_plots(self):
        """
        Generate all visualization plots.
        """
        print("\n" + "="*70)
        print(" GENERATING PHASE 2 VISUALIZATIONS")
        print("="*70 + "\n")
        
        self.load_data_and_models()
        
        self.plot_model_comparison()
        self.plot_predicted_vs_actual()
        self.plot_residual_analysis()
        self.plot_feature_importance()
        self.plot_error_by_depth()
        self.plot_error_by_sector()
        self.plot_feature_correlation_heatmap()
        self.plot_learning_curves()
        
        print("\n" + "="*70)
        print(" ALL VISUALIZATIONS COMPLETE!")
        print("="*70)
        print(f"\n✓ 8 publication-quality figures saved to:")
        print(f"  {self.output_dir}/")
        print("\nGenerated plots:")
        print("  01_model_comparison.png          - R² and RMSE comparison")
        print("  02_predicted_vs_actual.png        - Scatter plots for both models")
        print("  03_residual_analysis.png          - Residual plots and distributions")
        print("  04_feature_importance.png         - Feature importance rankings")
        print("  05_error_by_depth.png             - Error analysis by depth")
        print("  06_error_by_sector.png            - Error analysis by ocean sector")
        print("  07_correlation_heatmap.png        - Feature correlation matrix")
        print("  08_prediction_distribution.png    - Distribution comparison")
        print("\nThese figures are ready for your conference paper! 📊\n")


def main():
    """
    Main execution function.
    """
    visualizer = Phase2Visualizer(
        data_path='data/processed/southern_ocean_training.csv',
        results_dir='results/phase2',
        models_dir='models/phase2',
        output_dir='results/phase2/figures'
    )
    
    visualizer.generate_all_plots()


if __name__ == "__main__":
    main()
