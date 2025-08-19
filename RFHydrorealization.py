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
##################################################################
# Mean vs Severity
sim_cols = [c for c in grf_df.columns if c.startswith("GRF_sim_")]
grf_df["GRF_MEAN"] = grf_df[sim_cols].mean(axis=1)
df["GRF_MEAN"] = grf_df["GRF_MEAN"].values
# Valid points only
mask = (df["Veg_Code"] != "Urban") & df["GRF_MEAN"].notna() & df["Severity"].notna()
x = df.loc[mask, "Severity"]
y = df.loc[mask, "GRF_MEAN"]

# Correlation
pearson_r, _ = pearsonr(x, y)
spearman_r, _ = spearmanr(x, y)

# Plot
plt.figure(figsize=(7, 5))
plt.scatter(x, y, alpha=0.3, s=3, color="mediumseagreen")
plt.xlabel("Burn Severity")
plt.ylabel("Mean GRF Multiplier")
plt.title(f"GRF Mean vs Burn Severity\nPearson r = {pearson_r:.3f}, Spearman ρ = {spearman_r:.3f}")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.show()
##################################################################
# One realization vs Severity
realization = grf_samples[0, :]
mask = (df["Veg_Code"] != "Urban") & (realization > 0) & (df["Severity"] > 0)
x = np.log(df.loc[mask, "Severity"])
y = np.log(realization[mask])

pearson_r, _  = pearsonr(x, y)
spearman_r, _ = spearmanr(x, y)

plt.figure(figsize=(7, 5))
plt.scatter(x, y, alpha=0.3, s=3)
plt.xlabel("log(Burn Severity)")
plt.ylabel("log(GRF Multiplier)")
plt.title(f"Log-Log: One Realization\nPearson r = {pearson_r:.3f}, Spearman ρ = {spearman_r:.3f}")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.show()
##################################################################
# Heatmap: burn severity
mask_urban = df["Veg_Code"] == "Urban"
mask_valid = ~mask_urban
plt.figure(figsize=(8, 6), dpi=100)
sc = plt.scatter(
    df.loc[mask_valid, "POINT_X"],
    df.loc[mask_valid, "POINT_Y"],
    c=df.loc[mask_valid, "Severity"],
    cmap="Oranges", s=1, rasterized=True
)
plt.scatter(
    df.loc[mask_urban, "POINT_X"],
    df.loc[mask_urban, "POINT_Y"],
    color="black", s=0.5, label="Urban", rasterized=True
)
plt.colorbar(sc, label="Burn Severity")
plt.xlabel("X")
plt.ylabel("Y")
plt.axis("equal")
plt.title("Burn Severity Heatmap")
plt.legend(loc="lower right", markerscale=4)
plt.tight_layout()
plt.show()

# Heatmap: GRF Mean
# 1. Compute the per‐point mean across your Monte Carlo runs
sim_cols = [c for c in grf_df.columns if c.startswith("GRF_sim_")]
grf_df["GRF_MEAN"] = grf_df[sim_cols].mean(axis=1)

# 2. Merge it back into your original df (assuming same row order)
df["GRF_MEAN"] = grf_df["GRF_MEAN"].values
plt.figure(figsize=(8, 6), dpi=100)
sc = plt.scatter(
    df.loc[mask_valid, "POINT_X"],
    df.loc[mask_valid, "POINT_Y"],
    c=df.loc[mask_valid, "GRF_MEAN"],
    cmap="PuBu", s=1, rasterized=True,
    # vmin=0.2, vmax=1.2  # optional fixed scale
)
plt.scatter(
    df.loc[mask_urban, "POINT_X"],
    df.loc[mask_urban, "POINT_Y"],
    color="black", s=0.5, label="Urban", rasterized=True
)
plt.colorbar(sc, label="GRF Mean Multiplier")
plt.xlabel("X")
plt.ylabel("Y")
plt.axis("equal")
plt.title("GRF Mean (t = 0)")
plt.legend(loc="lower right", markerscale=4)
plt.tight_layout()
plt.show()


# Heatmap: GRF Std
# 1) Identify your sim columns
sim_cols = [c for c in grf_df.columns if c.startswith("GRF_sim_")]

# 2) Compute per‐point standard deviation across the Monte Carlo runs
grf_df["GRF_STD"] = grf_df[sim_cols].std(axis=1)

# 3) Bring it back into your main df (assuming rows line up)
df["GRF_STD"] = grf_df["GRF_STD"].values
# Get min and max values
vmin = df.loc[mask_valid, "GRF_STD"].min()
vmax = df.loc[mask_valid, "GRF_STD"].max()

plt.figure(figsize=(8, 6), dpi=100)
sc = plt.scatter(
    df.loc[mask_valid, "POINT_X"],
    df.loc[mask_valid, "POINT_Y"],
    c=df.loc[mask_valid, "GRF_STD"],
    cmap="viridis", vmin=vmin, vmax=vmax, s=1, rasterized=True
)
plt.scatter(
    df.loc[mask_urban, "POINT_X"],
    df.loc[mask_urban, "POINT_Y"],
    color="black", s=0.5, label="Urban", rasterized=True
)

# Add colorbar and annotate min/max
cbar = plt.colorbar(sc)
cbar.set_label("GRF Std Dev (across simulations)")

# Add text at top and bottom
cbar.ax.text(1.05, 1.0, f"max = {vmax:.2f}", transform=cbar.ax.transAxes,
             ha='left', va='top', fontsize=12)
cbar.ax.text(1.05, 0.0, f"min = {vmin:.2f}", transform=cbar.ax.transAxes,
             ha='left', va='bottom', fontsize=12)

plt.xlabel("X")
plt.ylabel("Y")
plt.axis("equal")
plt.title("GRF Std Dev (t = 0)")
plt.legend(loc="lower right", markerscale=4)
plt.tight_layout()
plt.show()

###################### Spatial Correlation ############################
from skgstat import Variogram

# Choose one realization
grf_realization = grf_samples[0, :]

# Filter valid (non-urban, non-NaN)
mask_valid = (df["Veg_Code"] != "Urban") & np.isfinite(grf_realization)
df_valid = df.loc[mask_valid].copy()
realization_vals = grf_realization[mask_valid]
coords_sampled = df_valid[["POINT_X", "POINT_Y"]].values

# Compute variogram
V = Variogram(
    coords_sampled,
    realization_vals,
    normalize=False,
    n_lags=10,
    maxlag=300,
    model='exponential'
)

# Plot
fig = V.plot(show=False, hist=False)
fig.suptitle("Empirical Variogram of One GRF Realization", fontsize=13)
plt.tight_layout()
plt.show()