import xarray as xr
import numpy as np
import os
import sys
import yaml
import pandas as pd
import gc
from datetime import datetime, timedelta
from concurrent.futures import ProcessPoolExecutor, as_completed

# ==============================================================================
# CONFIGURATION LOADER
# ==============================================================================
def load_task_config(full_cfg, var_key, level):
    task = full_cfg['active_task']
    year = str(task['year'])
    try:
        mapping = full_cfg['var_mapping'][var_key]
        mpas_var = f"{mapping['MPAS_prefix']}_{level}hPa"
        an_var = mapping.get('analysis_name', mapping['MPAS_prefix'])
        era5_var = mapping['ERA5']
    except KeyError:
        sys.exit(f"ERROR: Variable '{var_key}' configuration missing in YAML.")
        
    y_data = full_cfg['experiments'][year]
    return {
        'year': year,
        'level': level,
        'var_label': var_key,
        'mpas_var': mpas_var,
        'an_var': an_var,
        'era5_var': era5_var,
        'target_p': float(level * 100),
        'invariant_file': full_cfg['paths']['invariant_file'],
        'clim_dir': full_cfg['paths']['clim_dir'],
        'output_dir': y_data['output_dir'],
        'gfs_an_dir': y_data['gfs_an_dir'],
        'experiments': y_data['configs'], 
        'regions': full_cfg['regions'],
        'num_workers': task['num_workers'],
        'use_gfs': task['use_gfs_reference'],
        'gfs_prefix': "x1.2621442.init",
        'lead_times': task['lead_times'],
        'eval_start': task['eval_start'],
        'eval_end': task['eval_end']
    }

# ==============================================================================
# GET ANALYSIS FIELD
# ==============================================================================
def get_analysis_field(ds_an, an_key, zg_centers):
    """ Derives physical variables from native MPAS-JEDI or GFS analysis fields. """
    R_d, kappa, p0 = 287.05, 0.2857, 100000.0
    p_total = ds_an['pressure'].values[0] if 'pressure' in ds_an else \
              (ds_an['pressure_base'].values[0] + ds_an['pressure_p'].values[0])

    if an_key == 'theta':
        field = ds_an['theta'].values[0] * (p_total / p0)**kappa
    elif an_key == 'height':
        field = zg_centers
    elif an_key in ds_an:
        field = ds_an[an_key].values[0]
    else:
        raise ValueError(f"Analysis variable '{an_key}' not found in file.")

    return field, p_total

# ==============================================================================
# WORKER FUNCTION
# ==============================================================================
def process_date(init_date, exp_name, masks, area, clim_cache, cfg, zg_centers):
    results = []
    exp_cfg = cfg['experiments'][exp_name]
    init_str = init_date.strftime("%Y%m%d%H")

    for lead in cfg['lead_times']:
        valid_date = init_date + timedelta(hours=lead)
        v_ts, v_h = valid_date.strftime("%Y-%m-%d_%H.%M.%S"), valid_date.strftime("%H")
        fc_file = f"{exp_cfg['fc_dir']}/{init_str}/diag.{v_ts}.nc"

        # Determine Analysis File Path
        if cfg['use_gfs']:
            an_file = f"{cfg['gfs_an_dir']}/{valid_date:%Y%m%d%H}/{cfg['gfs_prefix']}.{v_ts}.nc"
        else:
            an_root = exp_cfg.get('an_dir', cfg['gfs_an_dir'])
            an_file = f"{an_root}/{valid_date:%Y%m%d%H}/an/{exp_cfg.get('an_prefix', 'an')}.{v_ts}.nc"
            if not os.path.exists(an_file):
                an_file = f"{an_root}/{valid_date:%Y%m%d%H}/{exp_cfg.get('an_prefix', 'an')}.{v_ts}.nc"

        if os.path.exists(fc_file) and os.path.exists(an_file):
            try:
                with xr.open_dataset(fc_file) as ds_f, xr.open_dataset(an_file) as ds_an:

                    # Get Forecast and Analysis Variables
                    F = ds_f[cfg['mpas_var']].values[0]
                    A_raw, p = get_analysis_field(ds_an, cfg['an_var'], zg_centers)

                    # Vertical Interpolation
                    log_p, log_target = np.log(p), np.log(cfg['target_p'])
                    idx1 = np.argmax(p < cfg['target_p'], axis=1)
                    idx0 = np.clip(idx1 - 1, 0, p.shape[1] - 1)
                    ii = np.arange(p.shape[0])

                    lp0, lp1 = log_p[ii, idx0], log_p[ii, idx1]
                    v0 = A_raw[ii, idx0] if A_raw.ndim > 1 else A_raw[idx0]
                    v1 = A_raw[ii, idx1] if A_raw.ndim > 1 else A_raw[idx1]

                    with np.errstate(divide='ignore', invalid='ignore'):
                        A = np.where(np.abs(lp1 - lp0) > 1e-9, v0 + (log_target - lp0) * (v1 - v0) / (lp1 - lp0), v0)

                    # Anomaly Calculation (Forecast/Analysis minus Climatology)
                    C = clim_cache[v_h]
                    
                    print(f"Mean Climatology (C):    {np.nanmean(C):.2f}")
                    print(f"Mean Forecast (F):    {np.nanmean(F):.2f}")
                    print(f"Mean Analysis (A):    {np.nanmean(A):.2f}")
                    f_prime, a_prime = F - C, A - C

                    # Regional Aggregation
                    for reg_name, mask in masks.items():
                        valid_m = np.isfinite(f_prime[mask]) & np.isfinite(a_prime[mask])
                        if not np.any(valid_m): continue
                        f_v, a_v, w_v = f_prime[mask][valid_m], a_prime[mask][valid_m], area[mask][valid_m]
                        w_v /= np.sum(w_v) # Normalize weights
                        f_dev, a_dev = f_v - np.sum(f_v * w_v), a_v - np.sum(a_v * w_v)
                        num = np.sum(f_dev * a_dev * w_v)
                        den = np.sqrt(np.sum(f_dev**2 * w_v) * np.sum(a_dev**2 * w_v))
                        results.append({"Region": reg_name, "Lead_Time": lead, "exp": exp_name, 
                                        "Init_Date": init_date, "acc": num/den if den > 1e-12 else np.nan})
                del F, A_raw, A, f_prime, a_prime
            except Exception: continue
    gc.collect()
    return results

# ==============================================================================
# MAIN DRIVER
# ==============================================================================
if __name__ == "__main__":
    CONFIG_FILE = "config_acc.yaml"
    with open(CONFIG_FILE, 'r') as f:
        full_cfg = yaml.safe_load(f)

    # 1. Load Static and Invariant for Lat/Lon and zgrid
    ds_static = xr.open_dataset(full_cfg['paths']['static_file'], decode_times=False)
    area = ds_static.areaCell.values
    lats, lons = np.rad2deg(ds_static.latCell.values), np.rad2deg(ds_static.lonCell.values) % 360 

    with xr.open_dataset(full_cfg['paths']['invariant_file'], decode_times=False) as ds_inv:
        zg = ds_inv['zgrid'].values
        zg_centers = 0.5 * (zg[:-1] + zg[1:]) if zg.ndim == 1 else 0.5 * (zg[:, :-1] + zg[:, 1:])

    masks = {k: (lats >= v['lat_min']) & (lats <= v['lat_max']) & (lons >= v['lon_min']) & (lons <= v['lon_max']) 
             for k, v in full_cfg['regions'].items()}

    # 2. Variable and Level Loops
    for var_key in full_cfg['active_task']['variables']:
        for level in full_cfg['active_task']['levels']:
            cfg = load_task_config(full_cfg, var_key, level)
            os.makedirs(cfg['output_dir'], exist_ok=True)

            # 3. Cache Climatology
            print(f"[*] Caching Climatology for {var_key} at {level}hPa...")
            clim_cache = {}
            for h in ["00", "06", "12", "18"]:
                c_path = os.path.join(cfg['clim_dir'], f"{h}00utc/mean_remapped_3D.nc")
                with xr.open_dataset(c_path) as ds_c:
                    clim_cache[h] = ds_c[cfg['era5_var']].sel(level=level, method="nearest").values[0]

            # 4. Prepare Task Queue
            print(f"\n--- TASK INITIALIZATION: {var_key.upper()} {level}hPa ---")
            start_dt = datetime.strptime(cfg['eval_start'], "%Y-%m-%d")
            end_dt = datetime.strptime(cfg['eval_end'], "%Y-%m-%d")
            num_days = (end_dt - start_dt).days + 1
            
            print(f"  > Global Period: {start_dt.date()} to {end_dt.date()} ({num_days} days)")
            task_list = [(start_dt + timedelta(days=x), name) for name in cfg['experiments'].keys() for x in range(num_days)]

            print(f"--- TOTAL RUNS TO PROCESS: {len(task_list)} ---\n")
            # 5. Parallel Execution
            raw_data = []
            with ProcessPoolExecutor(max_workers=cfg['num_workers']) as executor:
                futures = [executor.submit(process_date, d, e, masks, area, clim_cache, cfg, zg_centers) for d, e in task_list]
                for i, f in enumerate(as_completed(futures), 1):
                    raw_data.extend(f.result())
                    if i % 20 == 0: print(f"  > [{var_key}@{level}] Progress: {i}/{len(task_list)} runs completed")

            # 6. Aggregate and Export Statistics (Pivoted Wide Format)
            if raw_data:
                df = pd.DataFrame(raw_data)
                stats = df.groupby(['Region', 'Lead_Time', 'exp'])['acc'].agg(['mean', 'std', 'count']).reset_index()
                stats['ci95'] = 1.96 * (stats['std'] / np.sqrt(stats['count']))
                
                # Pivot Mean values: Region/Lead_Time as rows, Experiments as columns
                mean_pivot = stats.pivot(index=['Region', 'Lead_Time'], columns='exp', values='mean')

                # Pivot CI values with prefix
                ci_pivot = stats.pivot(index=['Region', 'Lead_Time'], columns='exp', values='ci95').add_prefix('ci_')
                df_final = pd.concat([mean_pivot, ci_pivot], axis=1).reset_index()

                out_file = f"ACC_{var_key}_{level}_{cfg['year']}.csv"
                df_final.to_csv(os.path.join(cfg['output_dir'], out_file), index=False)
                print(f"[SUCCESS] Written: {out_file}\n")

    print("\n[FINISH] All requested tasks completed.")
