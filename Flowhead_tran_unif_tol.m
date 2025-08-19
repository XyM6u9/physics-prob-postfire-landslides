function h = Flowhead_tran_unif_tol(L, L_1, gamaW, EMd, EMs, mu, theta_r, theta_s, k_s, alpha, m, htop, z, t, S, beta, analysmode, tolerance)
% ------------------------------------------------------------------------
% Purpose: Compute transient hydraulic head "h" at depth z and time t
%          under uniform root-zone sink and surface flux conditions.
% Inputs:
%   L          - Total soil column length (m)
%   L_1        - Depth of root-free zone (m)
%   gamaW      - Unit weight of water (kN/m^3)
%   EMd, EMs   - Elastic moduli (dry and saturated) used in deformation
%   mu         - Poisson's ratio
%   theta_r,s  - Residual and saturated water contents
%   k_s        - Saturated hydraulic conductivity (m/day)
%   alpha, m   - Gardner SWCC parameters (exponential model)
%   htop       - Applied head at the surface boundary (m)
%   z          - Depth at which head is evaluated (m)
%   t          - Elapsed time since initial condition (days)
%   S          - Uniform sink rate (root uptake per depth) (1/day)
%   beta       - Slope angle in radians
%   analysmode - 'coupled' or 'uncoupled' coupling with deformation
%   tolerance  - Convergence tolerance for fixed-point iteration
% Outputs:
%   h          - Computed hydraulic head (m)
% ------------------------------------------------------------------------
% analysmode='coupled';
% 1) Input validation: ensure positive k_s and alpha
assert(k_s > 0 && alpha > 0, 'k_s and alpha must be positive');

% 2) Precompute geometric and steady-state components:
L_2 = L - L_1;  % Root-zone thickness

% Steady-state solution h_ss(z):
c2 = (exp(alpha*htop) ...
    - exp(-alpha*L*cos(beta)) ...
    - S*L_2/(k_s*cos(beta)) ...
    - exp(-alpha*L_2*cos(beta))*S/(alpha*k_s*cos(beta)^2) ...
    + S/(alpha*k_s*cos(beta)^2)) ...
    / (1 - exp(-alpha*L*cos(beta)));
if z >= L_1
    % In the root-active zone:
    hss_bar = c2*(1 - exp(-alpha*z*cos(beta))) ...
        + S*(z - L_1)/(k_s*cos(beta)) ...
        + exp(-alpha*(z - L_1)*cos(beta)) * S/(alpha*k_s*cos(beta)^2) ...
        - S/(alpha*k_s*cos(beta)^2);
else
    % In the upper root-free zone:
    hss_bar = c2*(1 - exp(-alpha*z*cos(beta))); 
end

% 3) Transient coefficients (Fourier series amplitude weights):
gamma_1 = (exp(alpha*htop) ...
    - exp(-alpha*L*cos(beta)) ...
    - S*L_2/(k_s*cos(beta)) ...
    + S/(alpha*k_s*cos(beta)^2) ...
    - S*exp(-alpha*L_2*cos(beta))/(alpha*k_s*cos(beta)^2)) ...
    / (1 - exp(-alpha*L*cos(beta)));
gamma_2 = S/(k_s*cos(beta));
gamma_3 = S/(alpha*k_s*cos(beta)^2);

% 4) Iteratively solve coupled head-deformation fixed-point if required
if strcmp(analysmode, 'coupled')
    % Fixed-point iteration settings
    MAX_ITERS   = 5000;
    tol_rel     = 1e-4;
    tol_abs     = 1e-6;
    relax_ratio = 0.5;  % Under-relaxation factor
    h_old       = 0;

    for iter = 1:MAX_ITERS
        % 4a) Compute saturation Se from SWRC at previous head
        Se = SWRC(alpha, h_old);
        Se = min(max(Se, eps), 1 - eps);  % keep Se in (0,1)

        % 4b) Update effective stiffness c from elastic model
        E  = Elastic(EMd, EMs, Se, m);
        c  = (1/k_s) * (Se * gamaW * (1 + mu)*(1 - 2*mu)/(E*(1 - mu)) ...
                       + alpha * (theta_s - theta_r));

        % 4c) Compute transient Fourier sums for given c
        [sum_1,sum_2,sum_3,sum_4,sum_5,sum_6,sum_7,sum_8,sum_9,sum_10,sum_11,sum_12] = Sum_uniform_2(L, L_1, alpha, z, c, beta, t);

        % 4d) Assemble transient head correction h_hat(z,t)
        %    (combining gamma1, gamma2, gamma3 and sums(1:12))
        h_hat = assemble_h_hat(alpha, z, L, L_1, beta, c, gamma_1, gamma_2, gamma_3, sum_1,sum_2,sum_3,sum_4,sum_5,sum_6,sum_7,sum_8,sum_9,sum_10,sum_11,sum_12);
        
        % 4e) Fixed-point update: h_new = log transform
        arg = h_hat + hss_bar + exp(-alpha*cos(beta)*z);
        arg = max(arg, eps);
        assert(arg > 0, 'Negative argument in log');
        h_new = (1/alpha) * log(arg);

        % 4f) Under-relaxation and convergence check
        h = relax_ratio*h_new + (1-relax_ratio)*h_old;
        if abs(h - h_old) < max(tol_abs, tol_rel * abs(h))
            break;
        end
        h_old = h;
    end
    if iter == MAX_ITERS
        warning('Max iterations reached without convergence');
    end
else
    % Uncoupled mode: skip fixed-point and use linear c
    c    = alpha/k_s * (theta_s - theta_r);
    [sum_1,sum_2,sum_3,sum_4,sum_5,sum_6,sum_7,sum_8,sum_9,sum_10,sum_11,sum_12] = Sum_uniform_2(L, L_1, alpha, z, c, beta, t);
    h_hat = assemble_h_hat(alpha, z, L, L_1, beta, c, gamma_1, gamma_2, gamma_3, sum_1,sum_2,sum_3,sum_4,sum_5,sum_6,sum_7,sum_8,sum_9,sum_10,sum_11,sum_12);
    arg   = h_hat + hss_bar + exp(-alpha*cos(beta)*z);
    h     = (1/alpha) * log(arg);
end
end
