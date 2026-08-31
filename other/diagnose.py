import numpy as np
import scipy.io as sio

mat = sio.loadmat('data/raw/GLODAPv2.2023_Merged_Master_File.mat', squeeze_me=True)
AOU = mat['G2aou']
Region = mat['G2region']

for b_id, name in {1:"Atl", 8:"Pac", 16:"Ind", 4:"Arc"}.items():
    basin_aou = AOU[(Region == b_id) & (AOU > -900)]
    print(f"{name} AOU: Mean={np.mean(basin_aou):.3f}, Max={np.max(basin_aou):.3f}, 99th Pctl={np.percentile(basin_aou, 99):.3f}")