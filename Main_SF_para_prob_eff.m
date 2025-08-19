% --- Data Loading and Preallocation --------------------------------------
GridInfo = readtable('RasterT_Palisad4_SpatialJoin6_TableToExcel.xlsx', ...
    'Sheet', 'RasterT_Palisad4_SpatialJoin6');  % Load spatial join results

% Load GRF simulation tables from Excel
GRF_table  = readtable('GRF_Simulation_t32.xlsx');      % Each column: GRF_sim_1 ... GRF_sim_N
Root_table = readtable('Root_Simulation_t32.xlsx');     % Each column: Root_Cohesion_kPa_t0_run1 ... runN

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
slopeRad = deg2rad(slopeDeg);        % Slope angle (radians)

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

% --- Hydrological parameters and root-sink setup ----------------------------
L1        = 0.5;                                 % Depth up to no roots (m)
L2        = L - L1;
ET        = 1.2/24;                              % Evapotranspiration (mm/h)
FPAR      = 0.32;                                % Fraction for transpiration
Transp_val= ET * FPAR / 1000;                    % Convert mm to m

%%
% t         = 3;                                   % Time steps
% q_vec     = (1.57 * 0.0254) * ones(t,1);         % Infiltration rate % vector (m/h) %%benchmark
% q_vec     = (1.42 * 0.0254) * ones(t,1);         % Infiltration rate vector (m/h)

t         = 2;                                   % Time steps
% q_vec     = (1.89 * 0.0254) * ones(t,1);         % Infiltration rate vector (m/h)
q_vec     = (1.71 * 0.0254) * ones(t,1);         % Infiltration rate vector (m/h)

% t         = 6;                                   % Time steps
% q_vec     = (1.14 * 0.0254) * ones(t,1);         % Infiltration rate vector (m/h)
% q_vec     = (1.03 * 0.0254) * ones(t,1);         % Infiltration rate vector (m/h)

%%
Sink_vec  = Transp_val / L2 * ones(t,1);         % Sink term vector (m/h)
tolerance = 1e-4;                               % Solver tolerance

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
parpool('local_20core', 27);
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
        % Compute scaled conductivity m/h
        k_s = k_s_base(i) * grfRow(i) / 24;
        % Solve for pressure head profile
        try
            % Pressure head solver: use proper time-series vectors for q and Sink
            hM = case_1_hM_try_F_tol(L, L1, 2.00, gammaW, ...
                 EMd_arr(i), EMs_arr(i), 0.4, theta_r(i), theta_s(i), k_s, ...
                 alpha_arr(i), m_arr(i), slopeDeg(i), ...               % slope angle input
                 t, q_vec, Sink_vec, tolerance);
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
% Plotting code omitted for clarity

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
scatter(x, y, 3, values, 'filled',"square"); % Adjust marker size (20) as needed
cmap0 = [0.9 0.9 0.9];  
% cmap1 = hot(400);%12
% % cmap2 = flipud(bone(160));
% cmap2 = sky(3600);%40
cmap1 = hot(1200);%12
% % cmap2 = flipud(bone(160));
cmap2 = sky(8800);%40
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
percentage = sum(nonZeroValues < 1.2) / numel(nonZeroValues) * 100
%
figure(2)
% Extract coordinates and compute standard deviation
Coord = GridInfo{:, 6:7};  % Use curly braces to extract the actual data from table variables
compiled_std = std(All_SF1)'; % Assuming All_SF1 is a matrix where std is applied across the correct dimension
map_SF = [Coord compiled_std];
% Extract x, y, and scalar values
x = map_SF(:, 1);
y = map_SF(:, 2);
values = map_SF(:, 3);

% Set colormap: black for 0, sky for >0
numSkyColors = 400000;
cmap0 = [0.9 0.9 0.9];              % black for zero
cmap1 = hot(numSkyColors);   % sky colormap for positive values
cmap = [cmap0; cmap1];
colormap(cmap)

% Normalize color mapping: values from 0 to max expected (e.g., 4)
scatter(x, y, 14, values, 'filled', 's');
colorbar;
clim([prctile(compiled_std, 5) prctile(compiled_std, 95)]);  % sets clim
% clim([0 1]);  % sets clim
clim([0 0.2]);  % sets clim
xlabel('Northing (m)');
ylabel('Easting (m)');
box on
axis equal
nonZeroValues = values(values ~= 0);                      % Exclude zero elements
prc5=prctile(nonZeroValues, 5)
prc95=prctile(nonZeroValues, 50)
prc95=prctile(nonZeroValues, 95)