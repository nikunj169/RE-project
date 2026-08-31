import scipy.io as sio
import numpy as np

def preprocess_glodap(input_path, output_path):
    print(f"Loading {input_path}...")
    mat = sio.loadmat(input_path, squeeze_me=True)
    
    # Extract columns - Update keys based on your specific .mat structure
    S = mat['G2salinity']
    T = mat['G2temperature']
    AOU = mat['G2aou']
    # If you are validating against TCO2 or similar, extract it here:
    Y = mat.get('G2tco2', T) # Defaults to T if target not found
    
    # 1. Remove GLODAP null values (-999)
    mask = (S > 0) & (T > -10) & (AOU > -50) & (Y > -999)
    
    # 2. Critical Step: Avoid the Singularity
    # Your equation breaks at AOU = 3.054312. We must mask out values 
    # exactly equal to this to avoid 'inf'
    mask = mask & (np.abs(AOU - 3.054312) > 1e-4)
    
    S_clean = S[mask]
    T_clean = T[mask]
    AOU_clean = AOU[mask]
    Y_clean = Y[mask]

    # Save as compressed numpy file for speed
    np.savez(output_path, S=S_clean, T=T_clean, AOU=AOU_clean, Y=Y_clean)
    print(f"Processed data saved to {output_path}")
    print(f"Remaining data points: {len(S_clean)}")

if __name__ == "__main__":
    preprocess_glodap('data/raw/GLODAPv2.2023_Merged_Master_File.mat', 'processed_ocean_data.npz')