import xarray as xr
import numpy as np
import os
from dataclasses import dataclass
from typing import List
from datetime import datetime

@dataclass
class VariableName:
    era5_name: str
    mpas_name: str

BASE_VARS: List[VariableName] = [
    VariableName("Z", "Z"),
    VariableName("T", "T"),
    VariableName("U", "U"),
    VariableName("V", "V"),
    VariableName("SP", "SP"),
    VariableName("R", "R"),
    VariableName("Q", "Q"),
]

# --- Paths ---
STATIC_FILE = "/glade/campaign/mmm/parc/taosun/pandac/MPAS_GRAPH/x1.2621442.static.nc"
CLIM_BASE_PATH = "/glade/derecho/scratch/schwartz/MPAS/ERA5_forACC"
OUTPUT_BASE_DIR = "/glade/derecho/scratch/ivette/ERA5_remappedACC_MPAS"

G = 9.80665 

def get_climatology_path(valid_time: datetime) -> str:
    hour_str = valid_time.strftime("%H00utc")
    return os.path.join(CLIM_BASE_PATH, hour_str, "average", "mean.nc")

# 1. Load MPAS Mesh coordinates
ds_static = xr.open_dataset(STATIC_FILE, decode_times=False)
m_lat = np.rad2deg(ds_static.latCell).astype(np.float32)
m_lon = np.rad2deg(ds_static.lonCell).astype(np.float32) % 360

for h in [0, 6, 12, 18]:
    dummy_time = datetime(2024, 1, 1, h, 0) 
    cl_file = get_climatology_path(dummy_time)
    
    hour_folder = dummy_time.strftime("%H00utc")
    out_dir = os.path.join(OUTPUT_BASE_DIR, hour_folder)
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, "mean_remapped_3D.nc")
    
    if not os.path.exists(cl_file):
        print(f"Skipping: {cl_file} not found.")
        continue

    print(f"Processing {hour_folder}")
    with xr.open_dataset(cl_file, decode_times=False) as ds_cl:
        remapped_ds = xr.Dataset()
        
        # Copy original globals then update
        remapped_ds.attrs = ds_cl.attrs.copy()
        remapped_ds.attrs['grid_specification'] = f"Remapped to MPAS unstructured mesh ({len(m_lat)} cells)"
        remapped_ds.attrs['history'] = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}: Bilinear interpolation to MPAS mesh."

        if "level" in ds_cl.coords:
            remapped_ds.coords["level"] = ds_cl.coords["level"]
            remapped_ds.level.attrs = ds_cl.level.attrs.copy()
        
        encoding_dict = {}

        for var in BASE_VARS:
            if var.era5_name not in ds_cl:
                continue

            print(f"  Mapping {var.era5_name}...")
            
            data_interp = ds_cl[var.era5_name].interp(
                latitude=m_lat, 
                longitude=m_lon, 
                method="linear"
            ).compute().astype(np.float32)

            if var.era5_name == "Z":
                data_interp = data_interp / G
            
            remapped_ds[var.mpas_name] = data_interp
            remapped_ds[var.mpas_name].attrs = ds_cl[var.era5_name].attrs.copy()
            
            if var.era5_name == "Z":
                remapped_ds[var.mpas_name].attrs['units'] = 'm'
                remapped_ds[var.mpas_name].attrs['long_name'] = 'Geopotential Height'
                remapped_ds[var.mpas_name].attrs['standard_name'] = 'geopotential_height'

            orig_encoding = ds_cl[var.era5_name].encoding
            encoding_dict[var.mpas_name] = {
                'dtype': 'float32',
                '_FillValue': orig_encoding.get('_FillValue', None),
                'zlib': True,
                'complevel': 1
            }
            
        remapped_ds.to_netcdf(out_file, encoding=encoding_dict)
        print(f"DONE: {out_file}")
