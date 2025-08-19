import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt

# 1) Load full dataset
df = pd.read_excel("RasterT_Palisad4_SpatialJoin6_TableToExcel.xlsx")

# 2) Define vegetation classes to compute vs those to set NaN
valid_classes = {
    "Evergreen Needleleaf Forests",
    "Evergreen Broadleaf Forests",
    "Grasslands",
    "Shrublands",
    "Savannas"
}

# 3) Mark valid rows
is_valid = df["Veg_Code"].isin(valid_classes)
coords_valid = df.loc[is_valid, ["POINT_X", "POINT_Y"]].values
severity_valid = np.clip(df.loc[is_valid, "Severity"].values, 0, 1)
veg_valid = df.loc[is_valid, "Veg_Code"].values
n_valid = coords_valid.shape[0]
coords_urban = df.loc[df["Veg_Code"] == "Urban", ["POINT_X", "POINT_Y"]].values

# 4) Sidle (1991) parameters for each valid class
# base_params = {
#     'Evergreen Needleleaf Forests': {'k': 0.402, 'n': 0.647, 'a': 0.952, 'k1': 0.12, 'AC_inf': 12.5},
#     'Evergreen Broadleaf Forests':   {'k': 0.35,  'n': 0.65,  'a': 0.95,  'k1': 0.10, 'AC_inf': 10},
#     'Grasslands':                     {'k': 0.50,  'n': 0.60,  'a': 0.94,  'k1': 0.20, 'AC_inf': 3},
#     'Shrublands':                     {'k': 0.42,  'n': 0.64,  'a': 0.93,  'k1': 0.15, 'AC_inf': 6},
#     'Savannas':                       {'k': 0.45,  'n': 0.62,  'a': 0.92,  'k1': 0.13, 'AC_inf': 4.5},
# }
base_params = {
    'Evergreen Needleleaf Forests': {'k': 0.60, 'n': 0.85, 'a': 0.92, 'k1': 0.15, 'AC_inf': 12.5},
    'Evergreen Broadleaf Forests':   {'k': 0.38, 'n': 0.9, 'a': 0.95,  'k1': 0.15, 'AC_inf': 10},
    'Grasslands':                    {'k': 1.8,  'n': 1.2, 'a': 0.92,  'k1': 0.30, 'AC_inf': 3},
    'Shrublands':                    {'k': 0.7,  'n': 1.1, 'a': 0.93,  'k1': 0.18, 'AC_inf': 6},
    'Savannas':                      {'k': 0.6,  'n': 1.2, 'a': 0.91,  'k1': 0.20, 'AC_inf': 4.5},
}


# Compute b and c so that R(0)=0 and R(infty)=1
veg_params = {}
for name, p in base_params.items():
    a = p['a']
    b_new = a**2 / (1 - a)
    c_new = 1 - 1 / a
    veg_params[name] = {**p, 'b': b_new, 'c': c_new}

cv_map = {
    "Evergreen Needleleaf Forests": 0.60,
    "Evergreen Broadleaf Forests":  0.40,
    "Grasslands":                   0.70,
    "Shrublands":                   0.55,
    "Savannas":                     0.65
}

# 5) Precompute spatial neighbors & weights for valid points
nbrs = NearestNeighbors(n_neighbors=100).fit(coords_valid)
dists, idxs = nbrs.kneighbors(coords_valid)
# weights = np.exp(-dists / 500.0)
base_scale = 500
slope     = df["Slope"].values
mask_idx = np.where(is_valid)[0]
slope_valid = slope[mask_idx]  
scale_arr = base_scale * np.exp(-slope_valid / 30)
weights = np.exp(-dists / scale_arr[:, None]) 

# 6) Extract parameter arrays for valid points
k_arr = np.array([veg_params[v]["k"] for v in veg_valid])
n_arr = np.array([veg_params[v]["n"] for v in veg_valid])
a_arr = np.array([veg_params[v]["a"] for v in veg_valid])
b_arr = np.array([veg_params[v]["b"] for v in veg_valid])
c_arr = np.array([veg_params[v]["c"] for v in veg_valid])
k1_arr = np.array([veg_params[v]["k1"] for v in veg_valid])
AC_inf_arr = np.array([base_params[v]["AC_inf"] for v in veg_valid])
CV_arr = np.array([cv_map[v] for v in veg_valid])

# 7) Compute log-normal sampling parameters
sigma_ln_arr = np.sqrt(np.log(CV_arr**2 + 1))
mu_ln_arr = np.log(AC_inf_arr) - 0.5 * sigma_ln_arr**2

# 8) Time steps (0–12*80 months)
# t_months = np.arange(0,81)*12
# t_months = np.arange(0,13,3)*12
# t_months = np.arange(0,11,2)*12
# t_months = np.arange(0,121,30)
t_months = np.array([0, 8, 16, 32, 64, 128, 256, 512]) 
t_years = t_months / 12.0
n_times = len(t_months)

# 9) Allocate arrays for normalized infiltration and root cohesion
A_norm_valid = np.zeros((n_times, n_valid))
root_cohesion_valid = np.zeros((n_times, n_valid))

# 10) Compute for each time
for ti, t in enumerate(t_years):
    # compute normalized infiltration (no spatial smoothing on A)
    # C_a = np.clip(1 - severity_valid, 0.01, 0.99)
    C_a = np.clip(1 - severity_valid, 0.3, 0.99)
    D = C_a * np.exp(-k_arr * t**n_arr)
    R = 1 / (a_arr + b_arr * np.exp(-k1_arr * t)) + c_arr
    A_norm = D + R
    A_norm_valid[ti] = A_norm

    # sample raw C_max and then apply spatial smoothing
    C_max_raw = np.random.lognormal(mean=mu_ln_arr, sigma=sigma_ln_arr, size=n_valid)
    neigh_vals = C_max_raw[idxs]
    C_max_smooth = np.average(neigh_vals, axis=1, weights=weights)

    # compute root cohesion with spatially-correlated C_max
    root_cohesion_valid[ti] = A_norm * C_max_smooth

# 11) Build full-length outputs with NaN for invalid
n_total = df.shape[0]
# A_norm_full = np.full((n_times, n_total), np.nan)
# root_coh_full = np.full((n_times, n_total), np.nan)
A_norm_full    = np.zeros((n_times, n_total))
root_coh_full  = np.zeros((n_times, n_total))

mask_idx = np.where(is_valid)[0]
A_norm_full[:, mask_idx] = A_norm_valid
root_coh_full[:, mask_idx] = root_cohesion_valid

# # 12) Save t=0 snapshots back to DataFrame
# df["A_norm_t0"] = A_norm_full[0]
# df["Root_Cohesion_kPa_t0"] = root_coh_full[0]

######################################################
# 13) Monte Carlo sampling at t=0 for 100 runs
n_mc = 100  # number of Monte Carlo runs
A_norm_t0 = A_norm_full[0]
# root_coh_mc = np.full((n_mc, n_total), np.nan)
root_coh_mc  = np.zeros((n_mc, n_total)) 

for run in range(n_mc):
    raw = np.random.lognormal(mean=mu_ln_arr, sigma=sigma_ln_arr, size=n_valid)
    smooth = np.average(raw[idxs], axis=1, weights=weights)
    # tmp = np.full(n_total, np.nan)
    tmp     = np.zeros(n_total)
    # assign only valid points: multiply corresponding A_norm_t0 slice by smooth
    tmp_vals = A_norm_t0[mask_idx] * smooth
    tmp[mask_idx] = tmp_vals
    root_coh_mc[run] = tmp

# 14) Build output DataFrame and save to Excel
out_df = pd.DataFrame(index=df.index)
out_df['A_norm_t0'] = A_norm_t0
for run in range(n_mc):
    col_name = f"Root_Cohesion_kPa_t0_run{run+1}"
    out_df[col_name] = root_coh_mc[run]

filename = f"Root_Simulation_t0.xlsx"
out_df.to_excel(filename, index=False)
print(f"Saved simulation output to {filename}")
##############################################################################
# 15) Plot A_norm heatmaps at multiple times with Urban in black
# months_to_plot = [0, 12, 24, 60, 120, 240, 480, 960]
months_to_plot = [0, 8, 16, 32, 64, 128, 256, 512]
months_to_plot = list(t_months)
t_indices     = [np.where(t_months == m)[0][0] for m in months_to_plot]

# 计算色标归一化
vmin = np.nanmin(A_norm_valid[t_indices])
vmax = np.nanmax(A_norm_valid[t_indices])
norm = plt.Normalize(vmin=vmin, vmax=vmax)

# 创建子图
fig, axes = plt.subplots(4, 2, figsize=(10, 13), dpi=100,
                         gridspec_kw={"wspace": 0.05, "hspace": 0.15})
all_axes = axes.ravel().tolist()

# 科学计数格式化器
from matplotlib.ticker import ScalarFormatter, MaxNLocator
fmt = ScalarFormatter(useMathText=True)
fmt.set_powerlimits((0,0))

# 循环绘制
for ax, ti, m in zip(all_axes, t_indices, months_to_plot):
    # 主散点
    sc = ax.scatter(
        coords_valid[:, 0], coords_valid[:, 1],
        c=A_norm_valid[ti],
        cmap='viridis',
        norm=norm,
        s=5, marker='s', rasterized=True
    )
    # Urban 覆盖
    if coords_urban.size > 0:
        ax.scatter(
            coords_urban[:, 0], coords_urban[:, 1],
            color='k', s=5, marker='s', rasterized=True
        )

    # 设定坐标系、标签、formatter
    ax.set_aspect('equal', 'box')
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)
    ax.xaxis.set_major_locator(MaxNLocator(4))

    # 在内部左上角标注 t
    ax.text(
        0.05, 0.90,
        f"t = {m} months",
        transform=ax.transAxes,
        fontsize=12,
        va='top'
    )

# 只保留最外圈子图的刻度标签
for ax in all_axes:
    ax.label_outer()

# 添加单一 colorbar
cbar = fig.colorbar(
    sc, ax=all_axes,
    orientation='vertical',
    fraction=0.02, pad=0.04
)
cbar.set_label(r'$C^*$', fontsize=14)

plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()

################################
# —— Plot root cohesion heatmaps with Urban in black —— 
# specify the months you want to plot
# months_to_plot = [0, 12, 24, 60, 120, 240, 480, 960]
# find the corresponding time indices
t_indices = [np.where(t_months == m)[0][0] for m in months_to_plot]

# determine a common color range across all selected times
# 1) 计算色标归一化
vmin = np.nanmin(root_cohesion_valid[t_indices])
vmax = np.nanmax(root_cohesion_valid[t_indices])
norm = plt.Normalize(vmin=vmin, vmax=vmax)

# 2) 创建 2×4 的子图
fig, axes = plt.subplots(4, 2, figsize=(10, 13), dpi=100,
                         gridspec_kw={"wspace": 0.05, "hspace": 0.15})
all_axes = axes.ravel().tolist()

# 3) 科学计数格式化器
fmt = ScalarFormatter(useMathText=True)
fmt.set_powerlimits((0, 0))

# 4) 循环绘图
for ax, ti, m in zip(all_axes, t_indices, months_to_plot):
    # 主散点：root cohesion
    sc = ax.scatter(
        coords_valid[:, 0], coords_valid[:, 1],
        c=root_cohesion_valid[ti],
        cmap='PuBuGn',
        norm=norm,
        s=5, marker='s', rasterized=True
    )
    # Urban 覆盖
    if coords_urban.size > 0:
        ax.scatter(
            coords_urban[:, 0], coords_urban[:, 1],
            color='k', s=5, marker='s', rasterized=True
        )

    # 坐标系、刻度格式与刻度数量控制
    ax.set_aspect('equal', 'box')
    ax.set_xlabel('Easting (m)')
    ax.set_ylabel('Northing (m)')
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)
    ax.xaxis.set_major_locator(MaxNLocator(4))   # 最多 5 个刻度

    # 内部左上角标注时间
    ax.text(
        0.05, 0.90,
        f"t = {m} months",
        transform=ax.transAxes,
        fontsize=12,
        va='top'
    )

# 5) 只保留外圈的刻度标签
for ax in all_axes:
    ax.label_outer()

# 6) 添加共享 colorbar
cbar = fig.colorbar(
    sc, ax=all_axes,
    orientation='vertical',
    fraction=0.02, pad=0.04
)
cbar.set_label('Post-fire Root Cohesion (kPa)', fontsize=14)

# 7) 布局
plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()

##########
# ——— Scatter A_smooth vs. burn severity at t = 0 ———
ti0 = 0  # index for t = 0 months

# 1) raw normalized infil at t0 for valid points
A0 = A_norm_valid[ti0]

# 2) spatially smooth it using the same neighbors & weights
A0_smooth = np.average(A0[idxs], axis=1, weights=weights)

# 3) scatter
plt.figure(figsize=(8, 6))
plt.scatter(
    severity_valid,     # x = burn severity [0–1]
    A0_smooth,          # y = smoothed A_norm at t0
    s=20, alpha=0.7, edgecolor='k'
)
plt.xlabel('Burn Severity')
plt.ylabel('Spatially Smoothed $A_{norm}$ (t = 0)')
plt.title('Smoothed Normalized Infiltration vs. Burn Severity at t=0')
plt.grid(True)
plt.tight_layout()
plt.show()

#################
# — Compute mean and std across the 100 MC runs at t=0 — 
# mask_idx is the array of indices for valid (non-urban) points
# root_coh_mc has shape (n_mc, n_total)

# — Compute mean and std across the 100 MC runs at t=0 — 
mean_map = np.nanmean(root_coh_mc, axis=0)  # shape (n_total,)
std_map  = np.nanstd (root_coh_mc, axis=0)  # shape (n_total,)

# extract only the valid points
mask_idx    = np.where(is_valid)[0]
mean_valid  = mean_map[mask_idx]
std_valid   = std_map[mask_idx]

# — Figure 1: Mean root cohesion at t = 0 — 
plt.figure(figsize=(8, 6))
sc_mean = plt.scatter(
    coords_valid[:,0], coords_valid[:,1],
    c=mean_valid,
    cmap='viridis', s=5, marker='s'
)
plt.title('Mean Root Cohesion at t = 0')
plt.xlabel('POINT_X')
plt.ylabel('POINT_Y')
plt.gca().set_aspect('equal', 'box')
cbar = plt.colorbar(sc_mean, orientation='vertical', fraction=0.046, pad=0.04)
cbar.set_label('Mean Root Cohesion (kPa)')
plt.tight_layout()
plt.show()

# — Figure 2: Std dev of root cohesion at t = 0 — 
plt.figure(figsize=(8, 6))
sc_std = plt.scatter(
    coords_valid[:,0], coords_valid[:,1],
    c=std_valid,
    cmap='viridis', s=5, marker='s'
)
plt.title('Std Dev of Root Cohesion at t = 0')
plt.xlabel('POINT_X')
plt.ylabel('POINT_Y')
plt.gca().set_aspect('equal', 'box')
cbar = plt.colorbar(sc_std, orientation='vertical', fraction=0.046, pad=0.04)
cbar.set_label('Std Dev (kPa)')
plt.tight_layout()
plt.show()

###################### Spatial Correlation ############################
from skgstat import Variogram
import numpy as np
import matplotlib.pyplot as plt

# 1) take the full t=0 map (including NaN for Urban)
# real_full = root_cohesion_valid[0]   # shape = (n_total,)

# 2) mask out Urban & NaN via the original df
mask = is_valid 
coords_sampled = df.loc[mask, ["POINT_X", "POINT_Y"]].values
# values_sampled = real_full[mask]
values_sampled = root_cohesion_valid[0]

# 3) build & plot variogram
V = Variogram(
    coords_sampled,
    values_sampled,
    normalize=False,
    n_lags=5,
    maxlag=1000,
    model='exponential'
)

fig = V.plot(show=False, hist=False)
fig.suptitle("Empirical Variogram of $A_{norm}$ at t = 0", fontsize=13)
plt.tight_layout()
plt.show()
