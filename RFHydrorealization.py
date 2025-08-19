import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from scipy.stats import pearsonr, spearmanr

# Load and preprocess data
df = pd.read_excel("RasterT_Palisad4_SpatialJoin6_TableToExcel.xlsx")
df = df[["POINT_X", "POINT_Y", "Severity", "Slope", "Veg_Code", "Soil"]].dropna()
coords    = df[["POINT_X", "POINT_Y"]].values
severity  = df["Severity"].values
slope     = df["Slope"].values
veg_code  = df["Veg_Code"].values

# Simulation parameters
base_scale    = 500       # Spatial kernel base scale
k_neighbors   = 100       # Number of neighbors to consider
n_repeat      = 100       # Number of Monte Carlo runs
t_months      = 512#256         # Time since fire, in months
r_recovery    = 1.03      # Logistic recovery rate (per year)
tau_cv        = 5 * 12    # CV decay timescale (months)
mu_hydro      = -2.29     # log-mean for hydrophobic class
mu_philic     = 0.056     # log-mean for hydrophilic class
sigma_hydro   = 1.07      # log-std for hydrophobic class
sigma_philic  = 0.72      # log-std for hydrophilic class

# Build spatial neighbor model
nn_model = NearestNeighbors(n_neighbors=k_neighbors).fit(coords)
distances, indices = nn_model.kneighbors(coords)

# Prepare output array
n_points    = len(coords)
grf_samples = np.full((n_repeat, n_points), np.nan)

for run in range(n_repeat):
    for i in range(n_points):
        if veg_code[i] == "Urban":
            continue

        # Step 1: Get neighbor indices
        neighbors_idx = indices[i]
        
        # Step 2: Get neighbor attributes
        neighbor_severity = np.clip(severity[neighbors_idx], 0.01, 0.99)  # Ensure valid probability range
        neighbor_slope    = slope[neighbors_idx]
        local_dists       = distances[i]
        
        # Step 3: Compute spatial weights (distance-decay modulated by slope of cell i)
        weights = np.exp(-local_dists / (base_scale * np.exp(-slope[i] / 30)))

        # Step 4: Determine each neighbor's hydro status using its own severity
        p_hydro_neighbors = neighbor_severity
        is_hydro = np.random.rand(k_neighbors) < p_hydro_neighbors

        # Step 5: Assign lognormal params based on hydro status
        log_mu0  = np.where(is_hydro, mu_hydro, mu_philic)
        log_sig0 = np.where(is_hydro, sigma_hydro, sigma_philic)

        # Step 6: Recovery curve - normalized mean K*(t)
        K0      = np.exp(log_mu0)                             # initial K*
        a       = (1 - K0) / K0
        t_years = t_months / 12.0
        K_norm  = 1 / (1 + a * np.exp(-r_recovery * t_years))

        # Step 7: CV recovery
        cv0    = np.sqrt(np.exp(log_sig0**2) - 1)
        cv_inf = cv0 / 2.0
        cv_t   = cv_inf + (cv0 - cv_inf) * np.exp(-t_months / tau_cv)

        # Step 8: Invert CV to lognormal parameters at time t
        sigma_t = np.sqrt(np.log(cv_t**2 + 1))
        mu_t    = np.log(K_norm) - 0.5 * sigma_t**2

        # Step 9: Sample lognormal values for neighbors
        samples = np.random.lognormal(mean=mu_t, sigma=sigma_t)

        # Step 10: Weighted average + small noise
        value = np.average(samples, weights=weights)
        value += np.random.normal(scale=0.05)

        # Save result
        grf_samples[run, i] = value

        
# Create DataFrame from grf_samples (columns = simulation runs)
grf_df = pd.DataFrame(grf_samples.T, columns=[f"GRF_sim_{i+1}" for i in range(grf_samples.shape[0])])

# Concatenate coordinates
grf_df.insert(0, "POINT_Y", df["POINT_Y"].values)
grf_df.insert(0, "POINT_X", df["POINT_X"].values)

# Save to Excel
filename = f"GRF_Simulation_t{t_months}.xlsx"
grf_df.to_excel(filename, index=False)
print(f"Saved simulation output to {filename}")
