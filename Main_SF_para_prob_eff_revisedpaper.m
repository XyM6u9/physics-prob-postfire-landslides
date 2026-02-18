% --- Data Loading and Preallocation -------------------------------------- 
GridInfo = readtable('RasterT_Palisad4_SpatialJoin6_TableToExcel.xlsx', ...
    'Sheet', 'RasterT_Palisad4_SpatialJoin6');  % Load spatial join results

% Load GRF simulation tables from Excel
GRF_table  = readtable('GRF_Simulation_t8.xlsx');      % Each column: GRF_sim_1 ... GRF_sim_N
Root_table = readtable('Root_Simulation_t8.xlsx');     % Each column: Root_Cohesion_kPa_t0_run1 ... runN

% Sensitivy: Load Root simulation tables from Excel
% GRF_table  = readtable('GRF_Simulation_t16.xlsx');
% Root_table = readtable('Root_Sim_sensitivity_Evergreen_Forests_t16.xlsx'); %% %%%%%%%LINE 101%%%%%%%%%%%%%%
% Root_table = readtable('Root_Sim_sensitivity_Grasslands_t16.xlsx'); 
% Root_table = readtable('Root_Sim_sensitivity_ShrublandsSavannas_t16.xlsx');

% Convert category columns to string for grouping
GridInfo.Soil      = string(GridInfo.Soil);
GridInfo.LandCover = string(GridInfo.Veg_Code);

% Determine counts
[numCells, ~] = size(GridInfo);
numSim = 100;

% Preallocate result matrix (simulations x cells)
All_SF1 = zeros(numSim, numCells);

% --- Encode categories as numeric codes --------------------------------------
[soilCodes,  soilList]  = grp2idx(GridInfo.Soil);      % Soil type codes
[coverCodes, coverList] = grp2idx(GridInfo.LandCover); % Land cover codes
coverUrban = find(strcmp(coverList,'Urban'));
coverWater = find(strcmp(coverList,'Water'));

% Precompute slope angles in degrees and radians
slopeDeg = GridInfo.Slope;       % Slope angle (degrees)
slopeRad = deg2rad(slopeDeg);    % Slope angle (radians)

% --- Initialize soil property arrays ----------------------------------------
EMd_arr    = zeros(numCells,1);  % Dry elastic modulus kPa
EMs_arr    = zeros(numCells,1);  % Saturated elastic modulus kPa
theta_r    = zeros(numCells,1);  % Residual water content
theta_s    = zeros(numCells,1);  % Saturated water content
k_s_base   = zeros(numCells,1);  % Base saturated conductivity m/day
alpha_arr  = zeros(numCells,1);  % van Genuchten alpha 1/m
m_arr      = zeros(numCells,1);  % van Genuchten m
phi_deg    = zeros(numCells,1);  % Friction angle degrees
c_arr      = zeros(numCells,1);  % Cohesion kPa

for k = 1:numel(soilList)
    mask = soilCodes == k;
    switch soilList{k}
        case 'CL'
            EMd_arr(mask)=1870; EMs_arr(mask)=180;
            theta_r(mask)=0.06; theta_s(mask)=0.40;
            k_s_base(mask)=0.28; alpha_arr(mask)=0.90;
            m_arr(mask)=2.265; phi_deg(mask)=27; c_arr(mask)=1.0;
        case {'SC','SC-SM'}
            EMd_arr(mask)=3070; EMs_arr(mask)=675;
            theta_r(mask)=0.04; theta_s(mask)=0.39;
            k_s_base(mask)=0.55; alpha_arr(mask)=1.50;
            m_arr(mask)=1.140; phi_deg(mask)=32; c_arr(mask)=0;
        case {'SP','SP-SM'}
            EMd_arr(mask)=2980; EMs_arr(mask)=2883;
            theta_r(mask)=0.05; theta_s(mask)=0.375;
            k_s_base(mask)=24;   alpha_arr(mask)=2.60;
            m_arr(mask)=5.00;  phi_deg(mask)=37; c_arr(mask)=0;
        case 'SM'
            EMd_arr(mask)=3250; EMs_arr(mask)=3000;
            theta_r(mask)=0.04; theta_s(mask)=0.45;
            k_s_base(mask)=0.90; alpha_arr(mask)=0.20;
            m_arr(mask)=0.78;  phi_deg(mask)=32; c_arr(mask)=0.40;
        case 'GC'
            EMd_arr(mask)=4500; EMs_arr(mask)=1100;
            theta_r(mask)=0.01; theta_s(mask)=0.55;
            k_s_base(mask)=0.35; alpha_arr(mask)=0.10;
            m_arr(mask)=0.56;  phi_deg(mask)=30; c_arr(mask)=1.20;
    end
end

% Convert friction angles to radians
phi_arr = deg2rad(phi_deg);

gammaW = 9.8;   % Water unit weight kN/m^3
gammaS = 18;    % Soil unit weight kN/m^3

% --- Depth discretization ---------------------------------------------------
L = 1;                     % Profile depth (m)
Depth = (0:0.01:L)';       % Depth steps
nDepth = numel(Depth);

%%% NEW: SWR affected thickness (10 cm)
z_swr = 0.10;              % [m] 厚度
z_swr = min(z_swr, L);     % 防止比 profile 还厚

% --- Hydrological parameters and root-sink setup ----------------------------
L1        = 0.5;                                 % Depth up to no roots (m)
L2        = L - L1;
ET        = 1.2/24;                              % Evapotranspiration (mm/h)
FPAR      = 0.32;                                % Fraction for transpiration
Transp_val= ET * FPAR / 1000;                    % Convert mm to m

%%
t         = 3;                                   % Time steps
% q_vec     = (1.57 * 0.0254) * ones(t,1);         % Infiltration rate % vector (m/h) %%benchmark
q_vec     = (1.42 * 0.0254) * ones(t,1);         % Infiltration rate vector (m/h)

% t         = 2;                                   % Time steps
% q_vec     = (1.89 * 0.0254) * ones(t,1);         % Infiltration rate vector (m/h)
% q_vec     = (1.71 * 0.0254) * ones(t,1);         % Infiltration rate vector (m/h)

% t         = 6;                                   % Time steps
% q_vec     = (1.14 * 0.0254) * ones(t,1);         % Infiltration rate vector (m/h)
% q_vec     = (1.03 * 0.0254) * ones(t,1);         % Infiltration rate vector (m/h)

%% Sink term
Sink_vec  = Transp_val / L2 * ones(t,1);         % Sink term vector (m/h)
tolerance = 1e-4;                                % Solver tolerance

%%% NEW: TWI-based weighting for q (GridInfo 第 9 列)
TWI_raw = GridInfo{:,9};                         % 第九列是 TWI
isValidTWI = ~isnan(TWI_raw) & isfinite(TWI_raw);
TWI_valid  = TWI_raw(isValidTWI);

if ~isempty(TWI_valid)
    TWI_min = prctile(TWI_valid, 5);
    TWI_max = prctile(TWI_valid, 95);

    TWI_clip = TWI_raw;
    TWI_clip(TWI_clip < TWI_min) = TWI_min;
    TWI_clip(TWI_clip > TWI_max) = TWI_max;

    TWI_norm = (TWI_clip - TWI_min) ./ (TWI_max - TWI_min + eps);

    % 权重函数：高 TWI → 权重大 → q_eff 变大
    kappa = 1.0;   % 控制 TWI → q 的敏感度，可调参做敏感性
    w = exp(kappa * (TWI_norm - mean(TWI_norm(isValidTWI))));

    % 归一化，使得平均权重为 1（保证总水量守恒）
    w_eff = w / mean(w(isValidTWI));
else
    % 没有有效 TWI 时退化为均匀
    w_eff = ones(numCells,1);
end

% Preload GRF and Root matrices for faster access
GRF_mat  = zeros(numSim, numCells);
Root_mat = zeros(numSim, numCells);
for s = 1:numSim
    GRF_mat(s,:)  = GRF_table.(sprintf('GRF_sim_%d', s))';
    Root_mat(s,:) = Root_table.(sprintf('Root_Cohesion_kPa_t0_run%d', s))';
    
    % Sensitivity
    % GRF_mat(s,:)  = GRF_table.(sprintf('GRF_sim_%d', s))';
    % Root_mat(s,:) = Root_table.(sprintf('RootCoh_EvergreenForests_t16_run%d', s))';
    % Root_mat(s,:) = Root_table.(sprintf('RootCoh_Grasslands_t16_run%d', s))';
end

% Start parallel pool
parpool(24)%('local_20core', 27);
tic;
for s = 1:numSim
    % fprintf('Starting simulation %d of %d...\n', s, numSim);
    SF1 = zeros(1, numCells);
    grfRow  = GRF_mat(s,:);
    rootRow = Root_mat(s,:);
    parfor i = 1:numCells
        if mod(i,10000)==0
            fprintf("Simulation %d---%d/%d\n", s, i, numCells);
        end
        % Skip urban/water
        if coverCodes(i)==coverUrban || coverCodes(i)==coverWater
            SF1(i) = 0;
            continue;
        end

        % --- 1) 10 cm SWR 层 → 等效饱和导水率 Keff (m/h) -------------
        Ks0_day = k_s_base(i);         % [m/day] pre-fire Ks
        Ks0     = Ks0_day / 24;        % [m/h]

        Kstar       = grfRow(i);       % wildfire K* multiplier
        Kstar_safe  = max(Kstar, 1e-6);
        z1          = z_swr;           % SWR 层厚度
        z2          = L - z1;

        % 串联公式: 1/Keff = (z1/K1 + z2/K2)/L, K1=K*Ks0, K2=Ks0
        Keff_mult = L ./ ( z1 ./ Kstar_safe + z2 );   % 无量纲
        k_s       = Ks0 * Keff_mult;                  % [m/h] 等效 Ks

        % --- 2) TWI-based q 缩放 ----------------------------------------
        lambda_i = w_eff(i);                          % 该像元的 TWI 权重
        q_vec_i  = q_vec * lambda_i;                  % [t x 1] cell-specific q(t)

        % Solve for pressure head profile
        try
            % Pressure head solver: using q_vec_i and Sink_vec
            hM = case_1_hM_try_F_tol(L, L1, 2.00, gammaW, ...
                 EMd_arr(i), EMs_arr(i), 0.4, theta_r(i), theta_s(i), k_s, ...
                 alpha_arr(i), m_arr(i), slopeDeg(i), ...   % slope angle input
                 t, q_vec_i, Sink_vec, tolerance);
        catch
            hM = zeros(nDepth, 4);
        end

        % Vectorized FoS calculation
        Se1 = ones(size(hM));
        neg = hM<0;
        Se1(neg) = exp(alpha_arr(i) * hM(neg));
        rootRein = rootRow(i);
        denom = gammaS * sin(slopeRad(i)) * (L - Depth);
        SFmat = (c_arr(i) + rootRein - hM(:,2:end).*Se1(:,2:end).*tan(phi_arr(i))*gammaW) ...
                ./ denom + tan(phi_arr(i))/tan(slopeRad(i));
        vals = SFmat(SFmat>=0 & isfinite(SFmat));
        SF1(i) = min(vals);
    end
    All_SF1(s,:) = SF1;
end
toc;

%%
figure (1)
Coord = GridInfo{1:numCells, 6:7};
compiled_SF=mean(All_SF1)';
map_SF = [Coord compiled_SF];

% Assuming map_SF is an Nx3 matrix
x = map_SF(:, 1); % First column: X-coordinates
y = map_SF(:, 2); % Second column: Y-coordinates
values = map_SF(:, 3); % Third column: Values

% Create a scatter plot
scatter(x, y, 2, values, 'filled',"square"); % Adjust marker size (20) as needed
cmap0 = [0.9 0.9 0.9];  
% cmap1 = hot(400);%12
% % cmap2 = flipud(bone(160));
% cmap2 = sky(3600);%40
cmap1 = hot(1000);%12
% % cmap2 = flipud(bone(160));
cmap2 = sky(9000);%40
cmap = [cmap0;cmap1;cmap2];
% cmap = [cmap1;cmap2];
colormap(cmap)
clim([0 10]);
% clim([0 5]);
colorbar; % Add a colorbar to show the value scale
xlabel('Northing (m)');
ylabel('Easting (m)');
box on
axis equal; % Ensure equal scaling for x and y axes
nonZeroValues = values(values ~= 0);                      % Exclude zero elements
percentage = sum(nonZeroValues < 1.0) / numel(nonZeroValues) * 100
%


figure(2)
Coord = GridInfo{:, 6:7};
compiled_std = std(All_SF1)';
x = Coord(:,1);
y = Coord(:,2);
values = compiled_std;
% 分离 0 和 >0 的格子
zeroMask    = (values == 0);
nonZeroMask = ~zeroMask;
hold on
% 1) 先把 0 画成浅灰
scatter(x(zeroMask), y(zeroMask), 6, [0.85 0.85 0.85], 'filled', 's');

% 2) 再画非 0 的，用 colormap 映射
scatter(x(nonZeroMask), y(nonZeroMask), 6, values(nonZeroMask), ...
        'filled', 's');
% 用一个正常尺寸的 colormap
colormap(hot(256));         % 256 个颜色足够了
c = colorbar;
c.Label.String = 'Std(FS)';
% 颜色范围：你可以按数据设一个合理上限
v = values(nonZeroMask);
% clim([0, prctile(v, 95)]);  % 或者直接 [0 0.15] 之类
clim([0, 0.18]);
xlabel('Northing (m)');
ylabel('Easting (m)');
axis equal; box on;
nonZeroValues = values(values ~= 0);                      % Exclude zero elements
prc5=prctile(nonZeroValues, 5)
prc50=prctile(nonZeroValues, 50)
prc95=prctile(nonZeroValues, 95)

