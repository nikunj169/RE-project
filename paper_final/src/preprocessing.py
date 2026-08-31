import pandas as pd
import numpy as np
import scipy.io
import gsw
import os

class SouthernOceanDataPipeline:
    def __init__(self, filepath, save_path="data/processed/southern_ocean_training.csv"):
        self.filepath = filepath
        self.save_path = save_path
        self.df = None
        
        # Approximate Annual Atmospheric CO2 (Mauna Loa)
        self.co2_history = {
            1970: 325.68, 1975: 331.11, 1980: 338.75, 1985: 346.12,
            1990: 354.39, 1995: 360.82, 2000: 369.55, 2005: 379.80,
            2010: 389.90, 2015: 400.83, 2020: 414.24, 2023: 421.00
        }

    def load_mat_file(self):
        """Loads GLODAP .mat file and converts specific variables to DataFrame."""
        print(f"Loading GLODAP data from {self.filepath}...")
        
        try:
            mat = scipy.io.loadmat(self.filepath)
        except NotImplementedError:
            # Fallback for newer MATLAB v7.3 files (rare for GLODAP distribution but possible)
            import h5py 
            print("Detected v7.3 MAT file, switching to h5py...")
            mat = {}
            with h5py.File(self.filepath, 'r') as f:
                for k, v in f.items():
                    mat[k] = np.array(v)

        # Keys in GLODAP .mat usually match the CSV headers (e.g., 'G2latitude')
        # We need to extract them into a dictionary first
        
        # Define the variables we strictly need
        target_vars = [
            'G2latitude', 'G2longitude', 'G2year', 'G2month', 'G2depth', 'G2pressure',
            'G2temperature', 'G2salinity', 'G2salinityf',
            'G2oxygen', 'G2oxygenf',
            'G2fco2', 'G2fco2f',  # Fugacity
            'G2cruise', 'G2station'
        ]

        data_dict = {}
        
        # Check available keys (handling potential structural differences)
        available_keys = [k for k in mat.keys() if not k.startswith('__')]
        
        for var in target_vars:
            if var in available_keys:
                # Squeeze removes single-dimensional entries (e.g. (N,1) -> (N,))
                data_dict[var] = np.squeeze(mat[var])
            else:
                print(f"Warning: Variable {var} not found in .mat file.")
        
        self.df = pd.DataFrame(data_dict)
        print(f"Initial raw data size: {len(self.df)} rows")

    def filter_data(self):
        """Applies Southern Ocean + WOCE Flag filters."""
        if self.df is None: return

        # 1. Spatial Filter: Southern Ocean (Lat < -30)
        self.df = self.df[self.df['G2latitude'] <= -30].copy()
        print(f"Southern Ocean subset (Lat < -30): {len(self.df)} rows")

        # 2. Quality Control (WOCE Flags) - Keep Flag = 2 (Good)
        # Note: In .mat files, flags might be floats, ensure comparison works
        self.df = self.df[self.df['G2salinityf'].astype(int) == 2]
        self.df = self.df[self.df['G2oxygenf'].astype(int) == 2]
        
        # For fCO2, flag 2 is good. Sometimes 0 is used for calculated. 
        # We will stick to 2 for high quality training data.
        if 'G2fco2f' in self.df.columns:
            self.df = self.df[self.df['G2fco2f'].astype(int) == 2]
        
        # Drop NaNs in critical columns
        crit_cols = ['G2temperature', 'G2salinity', 'G2oxygen']
        if 'G2fco2' in self.df.columns: crit_cols.append('G2fco2')
        
        self.df.dropna(subset=crit_cols, inplace=True)
        print(f"High-quality data points remaining: {len(self.df)}")

    def feature_engineering(self):
        """Calculates Physics/Chemistry (Identical to previous steps)."""
        print("--- Starting Feature Engineering ---")

        # 1. Rename for clarity
        rename_map = {
            'G2latitude': 'latitude', 'G2longitude': 'longitude', 
            'G2year': 'year', 'G2depth': 'depth', 'G2pressure': 'pressure',
            'G2temperature': 'temperature', 'G2salinity': 'salinity',
            'G2oxygen': 'oxygen', 'G2fco2': 'fco2', 'G2cruise': 'cruise', 'G2station': 'station'
        }
        self.df.rename(columns=rename_map, inplace=True)

        # 2. TEOS-10 Variables (Using GSW)
        SA = gsw.SA_from_SP(self.df['salinity'].values, self.df['pressure'].values, 
                            self.df['longitude'].values, self.df['latitude'].values)
        CT = gsw.CT_from_t(SA, self.df['temperature'].values, self.df['pressure'].values)
        
        self.df['SA'] = SA
        self.df['CT'] = CT
        
        # 3. Calculate CO2(aq) from fCO2 (Weiss 1974)
        print("Calculating CO2(aq)...")
        Tk = self.df['temperature'] + 273.15
        S = self.df['salinity']
        
        # Coefficients (Weiss 1974)
        A1, A2, A3 = -60.2409, 93.4517, 23.3585
        B1, B2, B3 = 0.023517, -0.023656, 0.0047036
        
        ln_K0 = (A1 + A2 * (100/Tk) + A3 * np.log(Tk/100) + 
                 S * (B1 + B2 * (Tk/100) + B3 * (Tk/100)**2))
        
        K0 = np.exp(ln_K0) # mol/kg/atm
        self.df['co2_aq'] = self.df['fco2'] * K0 

        # 4. AOU (Water Mass Age Proxy)
        print("Calculating AOU...")
        O2_sat = gsw.O2sol(SA, CT, self.df['pressure'].values, 
                           self.df['longitude'].values, self.df['latitude'].values)
        self.df['aou'] = O2_sat - self.df['oxygen']

        # 5. Atmospheric Forcing
        years = np.array(list(self.co2_history.keys()))
        vals = np.array(list(self.co2_history.values()))
        self.df['c_atm'] = np.interp(self.df['year'], years, vals)

        # 6. Stratification (dS/dz)
        print("Calculating Stratification Gradients...")
        self.df.sort_values(by=['cruise', 'station', 'depth'], inplace=True)
        
        # Groupby is safer than shift logic when indexes are messy
        # Calculate diffs within groups
        g = self.df.groupby(['cruise', 'station'])
        d_sal = g['salinity'].diff()
        d_depth = g['depth'].diff()
        
        strat = d_sal / d_depth
        self.df['stratification_index'] = strat.fillna(0)
        
        # Clean up Infinite gradients (div by zero depth diffs)
        self.df['stratification_index'].replace([np.inf, -np.inf], 0, inplace=True)

        # 7. Sector Encoding
        lon = self.df['longitude']
        lon = np.where(lon > 180, lon - 360, lon)
        
        self.df['sector_atlantic'] = ((lon >= -70) & (lon < 20)).astype(int)
        self.df['sector_indian'] = ((lon >= 20) & (lon < 145)).astype(int)
        self.df['sector_pacific'] = 1 - (self.df['sector_atlantic'] + self.df['sector_indian'])

    def save_data(self):
        final_cols = [
            'year', 'latitude', 'longitude', 'depth', 'pressure',
            'temperature', 'salinity', 'aou', 'stratification_index',
            'c_atm', 'sector_atlantic', 'sector_indian', 'sector_pacific',
            'co2_aq'
        ]
        
        final_df = self.df[final_cols].replace([np.inf, -np.inf], np.nan).dropna()
        
        os.makedirs(os.path.dirname(self.save_path), exist_ok=True)
        final_df.to_csv(self.save_path, index=False)
        print(f"--- SUCCESS ---")
        print(f"Data saved to: {self.save_path}")
        print(f"Final Count: {len(final_df)} rows")

# --- Execution ---
if __name__ == "__main__":
    # Point to your .mat file
    MAT_PATH = "data/raw/GLODAPv2.2023_Merged_Master_File.mat" 
    
    if not os.path.exists(MAT_PATH):
        print(f"ERROR: File not found at {MAT_PATH}")
    else:
        pipeline = SouthernOceanDataPipeline(MAT_PATH)
        pipeline.load_mat_file()
        pipeline.filter_data()
        pipeline.feature_engineering()
        pipeline.save_data()