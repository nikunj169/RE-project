"""
Phase 3 - Equation Validation and Visualization
Uses the discovered symbolic equations to predict CO₂(aq) and compare with actual values

Equations:
- Full Model: (salinity/(3.054312 - aou) + temperature - (aou + 6.475268))**2 + 10.15233
- Simple Model: (11.368984 - temperature)*(aou - temperature) + 52.853294
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from mpl_toolkits.mplot3d import Axes3D
import os
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


class EquationValidator:
    """
    Validate and visualize the discovered symbolic equations.
    """
    
    def __init__(self, 
                 data_path='data/processed/southern_ocean_training.csv',
                 output_dir='results/phase3/validation'):
        """
        Initialize validator.
        """
        self.data_path = data_path
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        self.df = None
        self.scaler = None
        
    def load_data(self):
        """Load and prepare data."""
        print(f"Loading data from {self.data_path}...")
        self.df = pd.read_csv(self.data_path)
        print(f"✓ Loaded {len(self.df)} data points\n")
        
    def full_equation(self, salinity, temperature, aou, stratification_index, c_atm):
        """
        Full discovered equation (5 features).
        
        CO₂(aq) = (salinity/(3.054312 - aou) + temperature - (aou + 6.475268))**2 + 10.15233
        
        Note: Input features should be SCALED (standardized)
        """
        numerator = salinity
        denominator = 3.054312 - aou
        
        # Handle division by zero or near-zero
        denominator = np.where(np.abs(denominator) < 0.01, 0.01, denominator)
        
        term1 = numerator / denominator
        term2 = temperature
        term3 = aou + 6.475268
        
        inner = term1 + term2 - term3
        co2_aq = inner**2 + 10.15233
        
        return co2_aq
    
    def simple_equation(self, salinity, temperature, aou):
        """
        Simplified equation (3 features).
        
        CO₂(aq) = (11.368984 - temperature) × (aou - temperature) + 52.853294
        
        Note: Input features should be SCALED (standardized)
        """
        term1 = 11.368984 - temperature
        term2 = aou - temperature
        
        co2_aq = term1 * term2 + 52.853294
        
        return co2_aq
    
    def scale_features(self, df, feature_cols):
        """
        Scale features using StandardScaler.
        """
        X = df[feature_cols].values
        
        if self.scaler is None:
            self.scaler = StandardScaler()
            X_scaled = self.scaler.fit_transform(X)
        else:
            X_scaled = self.scaler.transform(X)
        
        return X_scaled
    
    def predict_full_model(self):
        """
        Use full equation to predict CO₂(aq) on entire dataset.
        """
        print("="*70)
        print(" FULL MODEL PREDICTION")
        print("="*70 + "\n")
        
        # Prepare features
        feature_cols = ['salinity', 'temperature', 'aou', 'stratification_index', 'c_atm']
        
        # Scale features
        X_scaled = self.scale_features(self.df, feature_cols)
        
        # Apply equation
        y_pred = self.full_equation(
            X_scaled[:, 0],  # salinity
            X_scaled[:, 1],  # temperature
            X_scaled[:, 2],  # aou
            X_scaled[:, 3],  # stratification_index
            X_scaled[:, 4]   # c_atm
        )
        
        # Actual values
        y_actual = self.df['co2_aq'].values
        
        # Metrics
        r2 = r2_score(y_actual, y_pred)
        rmse = np.sqrt(mean_squared_error(y_actual, y_pred))
        mae = mean_absolute_error(y_actual, y_pred)
        
        print(f"Full Model Performance:")
        print(f"  R² = {r2:.4f}")
        print(f"  RMSE = {rmse:.4f} μmol/kg")
        print(f"  MAE = {mae:.4f} μmol/kg\n")
        
        return y_pred, y_actual, r2, rmse
    
    def predict_simple_model(self):
        """
        Use simplified equation to predict CO₂(aq).
        """
        print("="*70)
        print(" SIMPLIFIED MODEL PREDICTION")
        print("="*70 + "\n")
        
        # Prepare features
        feature_cols = ['salinity', 'temperature', 'aou']
        
        # Reset scaler for simple model
        self.scaler = None
        X_scaled = self.scale_features(self.df, feature_cols)
        
        # Apply equation
        y_pred = self.simple_equation(
            X_scaled[:, 0],  # salinity
            X_scaled[:, 1],  # temperature
            X_scaled[:, 2]   # aou
        )
        
        # Actual values
        y_actual = self.df['co2_aq'].values
        
        # Metrics
        r2 = r2_score(y_actual, y_pred)
        rmse = np.sqrt(mean_squared_error(y_actual, y_pred))
        mae = mean_absolute_error(y_actual, y_pred)
        
        print(f"Simple Model Performance:")
        print(f"  R² = {r2:.4f}")
        print(f"  RMSE = {rmse:.4f} μmol/kg")
        print(f"  MAE = {mae:.4f} μmol/kg\n")
        
        return y_pred, y_actual, r2, rmse
    
    def plot_predicted_vs_actual(self, y_pred, y_actual, r2, model_name='Full Model'):
        """
        Create predicted vs actual scatter plot with non-linear trend analysis.
        """
        fig, axes = plt.subplots(1, 2, figsize=(18, 8))
        
        # LEFT PLOT: Traditional scatter
        ax1 = axes[0]
        
        # Scatter plot with density coloring
        scatter = ax1.scatter(y_actual, y_pred, 
                            alpha=0.4, s=15, c=y_actual, cmap='viridis', 
                            edgecolors='none')
        
        # Perfect prediction line
        min_val = min(y_actual.min(), y_pred.min())
        max_val = max(y_actual.max(), y_pred.max())
        ax1.plot([min_val, max_val], [min_val, max_val], 
               'r--', linewidth=3, label='Perfect Prediction', alpha=0.8)
        
        # Add polynomial trend line (degree 2 for quadratic)
        z = np.polyfit(y_actual, y_pred, 2)  # Degree 2 (quadratic)
        p = np.poly1d(z)
        x_trend = np.linspace(min_val, max_val, 300)
        ax1.plot(x_trend, p(x_trend), 
               'b-', linewidth=2.5, label=f'Quadratic Fit', alpha=0.7)
        
        ax1.set_xlabel('Actual CO₂(aq) (μmol/kg)', fontweight='bold', fontsize=13)
        ax1.set_ylabel('Predicted CO₂(aq) (μmol/kg)', fontweight='bold', fontsize=13)
        ax1.set_title(f'{model_name}: Predicted vs Actual\nR² = {r2:.4f}', 
                    fontweight='bold', fontsize=15)
        ax1.legend(loc='upper left', fontsize=11)
        ax1.grid(alpha=0.3)
        
        # Add statistics box
        rmse = np.sqrt(mean_squared_error(y_actual, y_pred))
        mae = mean_absolute_error(y_actual, y_pred)
        stats_text = f'R² = {r2:.4f}\nRMSE = {rmse:.2f}\nMAE = {mae:.2f}'
        ax1.text(0.05, 0.95, stats_text, transform=ax1.transAxes, 
               fontsize=11, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
        
        # Colorbar
        cbar = plt.colorbar(scatter, ax=ax1)
        cbar.set_label('Actual CO₂(aq)', rotation=270, labelpad=20, fontweight='bold')
        
        # RIGHT PLOT: Residual pattern showing non-linearity
        ax2 = axes[1]
        
        residuals = y_actual - y_pred
        
        # Hexbin plot to show density
        hexbin = ax2.hexbin(y_actual, residuals, gridsize=50, cmap='RdYlBu_r', 
                            mincnt=1, alpha=0.7)
        
        ax2.axhline(y=0, color='red', linestyle='--', linewidth=2.5, label='Zero Error')
        
        # Add polynomial trend to residuals (shows systematic bias)
        z_res = np.polyfit(y_actual, residuals, 2)
        p_res = np.poly1d(z_res)
        x_res = np.linspace(y_actual.min(), y_actual.max(), 300)
        ax2.plot(x_res, p_res(x_res), 
               'black', linewidth=2.5, label='Residual Trend (Quadratic)', alpha=0.8)
        
        ax2.set_xlabel('Actual CO₂(aq) (μmol/kg)', fontweight='bold', fontsize=13)
        ax2.set_ylabel('Residuals (Actual - Predicted)', fontweight='bold', fontsize=13)
        ax2.set_title(f'{model_name}: Residual Pattern\n(Shows Non-Linearity)', 
                    fontweight='bold', fontsize=15)
        ax2.legend(loc='upper left', fontsize=11)
        ax2.grid(alpha=0.3)
        
        # Colorbar for hexbin
        cbar2 = plt.colorbar(hexbin, ax=ax2)
        cbar2.set_label('Point Density', rotation=270, labelpad=20, fontweight='bold')
        
        plt.tight_layout()
        filename = f'01_{model_name.lower().replace(" ", "_")}_predicted_vs_actual.png'
        save_path = os.path.join(self.output_dir, filename)
        plt.savefig(save_path, bbox_inches='tight')
        print(f"✓ Saved: {save_path}")
        plt.close()
        
    def plot_residual_analysis(self, y_pred, y_actual, model_name='Full Model'):
        """
        Create residual analysis plots.
        """
        residuals = y_actual - y_pred
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # 1. Residuals vs Predicted
        ax1 = axes[0, 0]
        ax1.scatter(y_pred, residuals, alpha=0.3, s=10, color='steelblue')
        ax1.axhline(y=0, color='red', linestyle='--', linewidth=2)
        ax1.set_xlabel('Predicted CO₂(aq) (μmol/kg)', fontweight='bold')
        ax1.set_ylabel('Residuals (μmol/kg)', fontweight='bold')
        ax1.set_title('Residuals vs Predicted Values', fontweight='bold')
        ax1.grid(alpha=0.3)
        
        # 2. Residual histogram
        ax2 = axes[0, 1]
        ax2.hist(residuals, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        ax2.axvline(x=0, color='red', linestyle='--', linewidth=2)
        ax2.set_xlabel('Residuals (μmol/kg)', fontweight='bold')
        ax2.set_ylabel('Frequency', fontweight='bold')
        ax2.set_title('Residual Distribution', fontweight='bold')
        ax2.grid(alpha=0.3)
        
        # Add statistics
        mean_res = np.mean(residuals)
        std_res = np.std(residuals)
        ax2.text(0.05, 0.95, f'Mean: {mean_res:.3f}\nStd: {std_res:.3f}', 
                transform=ax2.transAxes, ha='left', va='top',
                fontsize=11, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
        
        # 3. Q-Q plot
        ax3 = axes[1, 0]
        from scipy import stats as scipy_stats
        scipy_stats.probplot(residuals, dist="norm", plot=ax3)
        ax3.set_title('Q-Q Plot (Normality Check)', fontweight='bold')
        ax3.grid(alpha=0.3)
        
        # 4. Absolute residuals vs predicted (heteroscedasticity check)
        ax4 = axes[1, 1]
        abs_residuals = np.abs(residuals)
        ax4.scatter(y_pred, abs_residuals, alpha=0.3, s=10, color='steelblue')
        
        # Add trend line
        z = np.polyfit(y_pred, abs_residuals, 1)
        p = np.poly1d(z)
        ax4.plot(sorted(y_pred), p(sorted(y_pred)), 'r-', linewidth=2, 
                label=f'Trend (slope={z[0]:.4f})')
        
        ax4.set_xlabel('Predicted CO₂(aq) (μmol/kg)', fontweight='bold')
        ax4.set_ylabel('|Residuals| (μmol/kg)', fontweight='bold')
        ax4.set_title('Absolute Residuals (Heteroscedasticity Check)', fontweight='bold')
        ax4.legend()
        ax4.grid(alpha=0.3)
        
        plt.suptitle(f'{model_name}: Residual Analysis', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        filename = f'02_{model_name.lower().replace(" ", "_")}_residual_analysis.png'
        save_path = os.path.join(self.output_dir, filename)
        plt.savefig(save_path, bbox_inches='tight')
        print(f"✓ Saved: {save_path}")
        plt.close()
        
    def plot_error_by_depth(self, y_pred, y_actual, model_name='Full Model'):
        """
        Analyze errors by depth bins.
        """
        depths = self.df['depth'].values
        errors = np.abs(y_actual - y_pred)
        
        # Define depth bins
        depth_bins = [0, 100, 500, 1000, 2000, 3000, 6000]
        depth_labels = ['0-100m', '100-500m', '500-1000m', 
                       '1000-2000m', '2000-3000m', '3000-6000m']
        
        depth_binned = pd.cut(depths, bins=depth_bins, labels=depth_labels)
        
        df_errors = pd.DataFrame({
            'depth_bin': depth_binned,
            'error': errors,
            'predicted': y_pred,
            'actual': y_actual
        })
        
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        # Box plot of errors
        ax1 = axes[0]
        df_errors.boxplot(column='error', by='depth_bin', ax=ax1, patch_artist=True)
        ax1.set_xlabel('Depth Range', fontweight='bold', fontsize=12)
        ax1.set_ylabel('Absolute Error (μmol/kg)', fontweight='bold', fontsize=12)
        ax1.set_title(f'{model_name}: Prediction Error by Depth', fontweight='bold', fontsize=13)
        plt.sca(ax1)
        plt.xticks(rotation=45, ha='right')
        ax1.get_figure().suptitle('')
        
        # Mean error by depth
        ax2 = axes[1]
        mean_errors = df_errors.groupby('depth_bin')['error'].mean()
        std_errors = df_errors.groupby('depth_bin')['error'].std()
        
        x_pos = np.arange(len(mean_errors))
        bars = ax2.bar(x_pos, mean_errors, yerr=std_errors, 
                      color='steelblue', alpha=0.7, capsize=5)
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(mean_errors.index, rotation=45, ha='right')
        ax2.set_xlabel('Depth Range', fontweight='bold', fontsize=12)
        ax2.set_ylabel('Mean Absolute Error (μmol/kg)', fontweight='bold', fontsize=12)
        ax2.set_title(f'{model_name}: Mean Error by Depth', fontweight='bold', fontsize=13)
        ax2.grid(axis='y', alpha=0.3)
        
        # Add value labels
        for i, bar in enumerate(bars):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.2f}',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        plt.tight_layout()
        filename = f'03_{model_name.lower().replace(" ", "_")}_error_by_depth.png'
        save_path = os.path.join(self.output_dir, filename)
        plt.savefig(save_path, bbox_inches='tight')
        print(f"✓ Saved: {save_path}")
        plt.close()
        
    def plot_error_by_sector(self, y_pred, y_actual, model_name='Full Model'):
        """
        Analyze errors by ocean sector.
        """
        # Determine sectors
        def get_sector(row):
            if row['sector_atlantic'] == 1:
                return 'Atlantic'
            elif row['sector_indian'] == 1:
                return 'Indian'
            else:
                return 'Pacific'
        
        sectors = self.df.apply(get_sector, axis=1)
        errors = np.abs(y_actual - y_pred)
        
        df_errors = pd.DataFrame({
            'sector': sectors,
            'error': errors,
            'predicted': y_pred,
            'actual': y_actual
        })
        
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        
        for idx, sector in enumerate(['Atlantic', 'Indian', 'Pacific']):
            ax = axes[idx]
            
            sector_data = df_errors[df_errors['sector'] == sector]
            
            # Scatter plot
            ax.scatter(sector_data['actual'], sector_data['predicted'], 
                      alpha=0.3, s=10, color='steelblue')
            
            # Perfect line
            min_val = sector_data['actual'].min()
            max_val = sector_data['actual'].max()
            ax.plot([min_val, max_val], [min_val, max_val], 
                   'r--', linewidth=2, label='Perfect Prediction')
            
            # Calculate R²
            r2_sector = r2_score(sector_data['actual'], sector_data['predicted'])
            rmse_sector = np.sqrt(mean_squared_error(sector_data['actual'], sector_data['predicted']))
            
            ax.set_xlabel('Actual CO₂(aq) (μmol/kg)', fontweight='bold')
            ax.set_ylabel('Predicted CO₂(aq) (μmol/kg)', fontweight='bold')
            ax.set_title(f'{sector} Sector\nR²={r2_sector:.4f}, RMSE={rmse_sector:.2f}', 
                        fontweight='bold')
            ax.legend(loc='upper left', fontsize=9)
            ax.grid(alpha=0.3)
            
            # Add sample size
            ax.text(0.05, 0.95, f'n={len(sector_data)}', 
                   transform=ax.transAxes, ha='left', va='top',
                   fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
        
        plt.suptitle(f'{model_name}: Performance by Ocean Sector', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        filename = f'04_{model_name.lower().replace(" ", "_")}_error_by_sector.png'
        save_path = os.path.join(self.output_dir, filename)
        plt.savefig(save_path, bbox_inches='tight')
        print(f"✓ Saved: {save_path}")
        plt.close()
        
    def plot_equation_surface(self, model_name='Full Model'):
        """
        Create 3D surface plots showing the non-linear equation behavior.
        """
        print(f"\nGenerating 3D surface plot for {model_name}...")
        
        if model_name == 'Full Model':
            # For full model: show Salinity vs AOU vs CO2
            feature_cols = ['salinity', 'temperature', 'aou', 'stratification_index', 'c_atm']
            self.scaler = None
            X_scaled = self.scale_features(self.df, feature_cols)
            
            # Create grid for salinity and AOU (scaled)
            sal_scaled = X_scaled[:, 0]
            aou_scaled = X_scaled[:, 2]
            temp_mean = X_scaled[:, 1].mean()
            strat_mean = X_scaled[:, 3].mean()
            catm_mean = X_scaled[:, 4].mean()
            
            # Create meshgrid
            sal_range = np.linspace(sal_scaled.min(), sal_scaled.max(), 50)
            aou_range = np.linspace(aou_scaled.min(), aou_scaled.max(), 50)
            SAL, AOU = np.meshgrid(sal_range, aou_range)
            
            # Calculate CO2 for each grid point
            CO2 = np.zeros_like(SAL)
            for i in range(SAL.shape[0]):
                for j in range(SAL.shape[1]):
                    CO2[i, j] = self.full_equation(
                        SAL[i, j], temp_mean, AOU[i, j], strat_mean, catm_mean
                    )
            
            xlabel = 'Salinity (scaled)'
            ylabel = 'AOU (scaled)'
            
        else:  # Simple model
            feature_cols = ['salinity', 'temperature', 'aou']
            self.scaler = None
            X_scaled = self.scale_features(self.df, feature_cols)
            
            # Create grid for temperature and AOU
            temp_scaled = X_scaled[:, 1]
            aou_scaled = X_scaled[:, 2]
            
            temp_range = np.linspace(temp_scaled.min(), temp_scaled.max(), 50)
            aou_range = np.linspace(aou_scaled.min(), aou_scaled.max(), 50)
            TEMP, AOU = np.meshgrid(temp_range, aou_range)
            
            # Calculate CO2
            CO2 = np.zeros_like(TEMP)
            for i in range(TEMP.shape[0]):
                for j in range(TEMP.shape[1]):
                    CO2[i, j] = self.simple_equation(
                        0, TEMP[i, j], AOU[i, j]  # salinity not used
                    )
            
            SAL = TEMP  # For consistent variable naming
            xlabel = 'Temperature (scaled)'
            ylabel = 'AOU (scaled)'
        
        # Create 3D plot
        fig = plt.figure(figsize=(16, 6))
        
        # Surface plot
        ax1 = fig.add_subplot(121, projection='3d')
        surf = ax1.plot_surface(SAL, AOU, CO2, cmap='viridis', 
                               alpha=0.8, edgecolor='none', antialiased=True)
        ax1.set_xlabel(xlabel, fontweight='bold', fontsize=11)
        ax1.set_ylabel(ylabel, fontweight='bold', fontsize=11)
        ax1.set_zlabel('CO₂(aq) (μmol/kg)', fontweight='bold', fontsize=11)
        ax1.set_title(f'{model_name}: 3D Surface\n(Non-Linear Parabolic Shape)', 
                     fontweight='bold', fontsize=13)
        ax1.view_init(elev=25, azim=45)
        fig.colorbar(surf, ax=ax1, shrink=0.5, aspect=5)
        
        # Contour plot
        ax2 = fig.add_subplot(122)
        contour = ax2.contourf(SAL, AOU, CO2, levels=20, cmap='viridis', alpha=0.8)
        contour_lines = ax2.contour(SAL, AOU, CO2, levels=10, colors='black', 
                                    linewidths=0.5, alpha=0.4)
        ax2.clabel(contour_lines, inline=True, fontsize=8)
        ax2.set_xlabel(xlabel, fontweight='bold', fontsize=11)
        ax2.set_ylabel(ylabel, fontweight='bold', fontsize=11)
        ax2.set_title(f'{model_name}: Contour Plot\n(Iso-CO₂ Lines)', 
                     fontweight='bold', fontsize=13)
        fig.colorbar(contour, ax=ax2, label='CO₂(aq) (μmol/kg)')
        
        plt.tight_layout()
        filename = f'06_{model_name.lower().replace(" ", "_")}_3d_surface.png'
        save_path = os.path.join(self.output_dir, filename)
        plt.savefig(save_path, bbox_inches='tight', dpi=200)
        print(f"✓ Saved: {save_path}")
        plt.close()
        
    def compare_both_models(self, y_pred_full, y_pred_simple, y_actual):
        """
        Compare full and simplified models side by side.
        """
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))
        
        models = [
            ('Full Model', y_pred_full, '#3498db'),
            ('Simple Model', y_pred_simple, '#e74c3c')
        ]
        
        for idx, (name, y_pred, color) in enumerate(models):
            ax = axes[idx]
            
            r2 = r2_score(y_actual, y_pred)
            rmse = np.sqrt(mean_squared_error(y_actual, y_pred))
            
            # Scatter
            ax.scatter(y_actual, y_pred, alpha=0.3, s=10, color=color)
            
            # Perfect line
            min_val = min(y_actual.min(), y_pred.min())
            max_val = max(y_actual.max(), y_pred.max())
            ax.plot([min_val, max_val], [min_val, max_val], 
                   'r--', linewidth=2, label='Perfect Prediction')
            
            # Quadratic fit
            z = np.polyfit(y_actual, y_pred, 2)
            p = np.poly1d(z)
            x_trend = np.linspace(min_val, max_val, 300)
            ax.plot(x_trend, p(x_trend), 
                   'black', linewidth=2, label='Quadratic Fit', alpha=0.6)
            
            ax.set_xlabel('Actual CO₂(aq) (μmol/kg)', fontweight='bold', fontsize=12)
            ax.set_ylabel('Predicted CO₂(aq) (μmol/kg)', fontweight='bold', fontsize=12)
            ax.set_title(f'{name}\nR²={r2:.4f}, RMSE={rmse:.2f}', 
                        fontweight='bold', fontsize=13)
            ax.legend(loc='upper left')
            ax.grid(alpha=0.3)
            
            # Stats box
            stats = f'R² = {r2:.4f}\nRMSE = {rmse:.2f}\nMAE = {mean_absolute_error(y_actual, y_pred):.2f}'
            ax.text(0.95, 0.05, stats, transform=ax.transAxes, 
                   ha='right', va='bottom', fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
        
        plt.suptitle('Comparison: Full vs Simplified Equation', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        save_path = os.path.join(self.output_dir, '05_model_comparison.png')
        plt.savefig(save_path, bbox_inches='tight')
        print(f"✓ Saved: {save_path}")
        plt.close()
        
    def run_complete_validation(self):
        """
        Run complete validation workflow.
        """
        print("\n" + "="*70)
        print(" EQUATION VALIDATION & VISUALIZATION")
        print(" Using Discovered Symbolic Equations")
        print("="*70)
        print("\nFull Equation:")
        print("  CO₂(aq) = (salinity/(3.054 - aou) + temp - (aou + 6.475))² + 10.152")
        print("\nSimple Equation:")
        print("  CO₂(aq) = (11.369 - temp) × (aou - temp) + 52.853")
        print("="*70 + "\n")
        
        self.load_data()
        
        # Full model
        print("\n" + "="*70)
        print(" VALIDATING FULL MODEL")
        print("="*70)
        y_pred_full, y_actual_full, r2_full, rmse_full = self.predict_full_model()
        
        self.plot_predicted_vs_actual(y_pred_full, y_actual_full, r2_full, 'Full Model')
        self.plot_residual_analysis(y_pred_full, y_actual_full, 'Full Model')
        self.plot_error_by_depth(y_pred_full, y_actual_full, 'Full Model')
        self.plot_error_by_sector(y_pred_full, y_actual_full, 'Full Model')
        self.plot_equation_surface('Full Model')
        
        # Simple model
        print("\n" + "="*70)
        print(" VALIDATING SIMPLIFIED MODEL")
        print("="*70)
        y_pred_simple, y_actual_simple, r2_simple, rmse_simple = self.predict_simple_model()
        
        self.plot_predicted_vs_actual(y_pred_simple, y_actual_simple, r2_simple, 'Simple Model')
        self.plot_residual_analysis(y_pred_simple, y_actual_simple, 'Simple Model')
        self.plot_error_by_depth(y_pred_simple, y_actual_simple, 'Simple Model')
        self.plot_error_by_sector(y_pred_simple, y_actual_simple, 'Simple Model')
        self.plot_equation_surface('Simple Model')
        
        # Comparison
        print("\n" + "="*70)
        print(" COMPARING BOTH MODELS")
        print("="*70 + "\n")
        self.compare_both_models(y_pred_full, y_pred_simple, y_actual_full)
        
        print("\n" + "="*70)
        print(" ✓ VALIDATION COMPLETE!")
        print("="*70)
        print(f"\nAll visualizations saved to: {self.output_dir}/")
        print("\n📊 Generated 11 Publication-Quality Figures:")
        print("\n  FULL MODEL:")
        print("    01_full_model_predicted_vs_actual.png     (with quadratic fit)")
        print("    02_full_model_residual_analysis.png       (4-panel diagnostics)")
        print("    03_full_model_error_by_depth.png          (depth stratification)")
        print("    04_full_model_error_by_sector.png         (Atlantic/Indian/Pacific)")
        print("    06_full_model_3d_surface.png              (parabolic surface)")
        print("\n  SIMPLE MODEL:")
        print("    01_simple_model_predicted_vs_actual.png")
        print("    02_simple_model_residual_analysis.png")
        print("    03_simple_model_error_by_depth.png")
        print("    04_simple_model_error_by_sector.png")
        print("    06_simple_model_3d_surface.png")
        print("\n  COMPARISON:")
        print("    05_model_comparison.png                   (side-by-side)")
        
        print("\n" + "="*70)
        print(" 🎯 KEY FINDINGS")
        print("="*70)
        print(f"\n  Full Model:   R² = {r2_full:.4f}, RMSE = {rmse_full:.2f} μmol/kg")
        print(f"  Simple Model: R² = {r2_simple:.4f}, RMSE = {rmse_simple:.2f} μmol/kg")
        print("\n  ✅ Non-linear (squared) relationships visualized")
        print("  ✅ Equations validated across depths and sectors")
        print("  ✅ Ready for conference paper!\n")


def main():
    """
    Main execution.
    """
    validator = EquationValidator(
        data_path='data/processed/southern_ocean_training.csv',
        output_dir='results/phase3/validation'
    )
    
    validator.run_complete_validation()


if __name__ == "__main__":
    main()
