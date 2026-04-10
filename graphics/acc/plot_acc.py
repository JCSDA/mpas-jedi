import yaml
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import numpy as np

def run_plotting_suite(config_path="config_acc.yaml"):
    # 1. Load YAML Configuration
    if not os.path.exists(config_path):
        print(f"ERROR: Config file '{config_path}' not found.")
        return

    with open(config_path, 'r') as f:
        full_cfg = yaml.safe_load(f)
    
    # 2. Extract Metadata and Paths
    task = full_cfg['active_task']
    year = str(task['year'])
    vars_to_plot = task['variables']
    levels_to_plot = task['levels']
    
    # Unit conversion for x-axis (Hours to Days)
    lead_times_days = np.array(task['lead_times']) / 24.0

    try:
        output_dir = full_cfg['experiments'][year]['output_dir']
    except KeyError:
        print(f"ERROR: Could not find output_dir for year {year} in YAML.")
        return

    print(f"[*] Processing plots in: {output_dir}")

    # 3. Iterate through Variables and Levels
    for var in vars_to_plot:
        for lev in levels_to_plot:
            csv_name = f"ACC_{var}_{lev}_{year}.csv"
            csv_path = os.path.join(output_dir, csv_name)
            
            if not os.path.exists(csv_path):
                print(f"  [!] Skipped: {csv_name} (Not found)")
                continue
                
            df = pd.read_csv(csv_path)

            # 4. Generate a plot for each Region
            for region in df['Region'].unique():
                reg_df = df[df['Region'] == region].copy()
                reg_df['Lead_Days'] = reg_df['Lead_Time'] / 24.0
                
                plt.figure(figsize=full_cfg.get('plotting', {}).get('fig_size', [10, 6]))
                
                plot_styles = full_cfg.get('plotting', {}).get('exp_styles', {})
                # Get experiment columns (exclude metadata and CI columns)
                exp_cols = [c for c in df.columns if c not in ['Region', 'Lead_Time'] and not c.startswith('ci_')]
                
                for exp in exp_cols:
                    # Get style from YAML; use key 'exp' as the label
                    style = plot_styles.get(exp, {})
                    color = style.get('color', None)
                    ls = style.get('ls', '-')
                    marker = style.get('marker', 'o')
                    
                    plt.plot(reg_df['Lead_Days'], reg_df[exp], 
                             label=exp,
                             color=color, 
                             linestyle=ls, 
                             marker=marker, 
                             markersize=5, 
                             linewidth=2.0)
                    
                    # Add Shaded Confidence Intervals
                    ci_col = f"ci_{exp}"
                    if ci_col in reg_df.columns:
                        plt.fill_between(reg_df['Lead_Days'], 
                                         reg_df[exp] - reg_df[ci_col], 
                                         reg_df[exp] + reg_df[ci_col], 
                                         color=color, alpha=0.15)

                # 5.Title and axes labels
                use_gfs = task['use_gfs_reference']
                ref = "GFS" if use_gfs else "MPAS-JEDI"
                plt.title(f"{var.upper()} @ {lev}hPa | {region} ({year}) | Analysis: {ref}", fontweight='bold', fontsize=14)
                plt.xlabel("Lead Time (days)", fontsize=12)
                plt.ylabel("Anomaly Correlation Coefficient", fontsize=12)
                
                plt.ylim(0.4, 1.02)
                plt.xticks(lead_times_days)
                plt.grid(True, linestyle='--', alpha=0.4)
                
                # Plot the 0.6 reference line (Standard for NWP skill)
                #plt.axhline(0.6, color='red', linestyle=':', alpha=0.5, label="Skill Limit (0.6)")
                
                plt.legend(loc='lower left', fontsize='medium', ncol=2, frameon=True)

                # 6. Save Plot
                save_name = f"ACC_{var}_{lev}_{region}_{year}.png"
                plt.savefig(os.path.join(output_dir, save_name), dpi=300, bbox_inches='tight')
                plt.close()
                print(f"  [+] Saved: {save_name}")

if __name__ == "__main__":
    run_plotting_suite()
