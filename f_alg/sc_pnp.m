function [X_est, time_est, S_est, C_est] = sc_pnp(X_observed, Omega, R, Denoiser, varargin)

% -------------------------------------------------------------------------
% --- Problem ---
% \minimize_{s_r \geq 0, c_r \geq 0}  
%       \|  O \oast (Y - \sum_{r=1}^R S_r \circ c_r)  \|_F^2
%       + \lambda \sum_{r=1}^R \phi(s_r) + \xi \sum_{r=1}^R c_r^T c_r
%
% with 
% O     -> measurement tensor (fiberwise), with O(i,j,:) = 1 being observed,
%          0 otherwise
% Y     -> observed Radio map
% S_r   -> r-th slf
% c_r   -> r-th psd
% \phi(s_r)     -> denoiser asscociated regularizer
% lambda, xi    -> regularization parameters
%                  lambda = rho * Noise_sigma, with Noise_sigma^2 denoting 
%                  the noise power used in denoisers, and rho denoting the
%                  augmented Lagrangian parameter
%
% --- Reference ---
% "Radio Map Estimation via Latent Domain Plug-and-Play Denoising",
% Le Xu, Lei Cheng, Junting Chen, Wenqiang Pu, and Xiao Fu
%
% Written by Le XU, updated by Nov 26, 2025
% -------------------------------------------------------------------------

%% Parse parameters
p = inputParser;
addOptional(p,'Denoiser',@isstring);

addOptional(p,'max_iter',30,@isnumeric);
addOptional(p,'threshold',1e-2,@isnumeric);
addOptional(p,'X_true',X_observed,@isnumeric);
addOptional(p,'S_true',0,@isnumeric); % for test only
addOptional(p,'C_true',0,@isnumeric); % for test only
addOptional(p,'BuildingMask',0,@isnumeric);  
addOptional(p,'NLM_MEMORY_EFF',false,@islogical)

switch Denoiser
    case "BM3D"
        defaultPar.xi = 1e-7;
        defaultPar.rho = 0.002;
        defaultPar.lambda = 20 * defaultPar.rho;
    case "DPIRNet"
        defaultPar.xi = 1e-7;
        defaultPar.rho = 0.01;
        defaultPar.lambda = 2.5e-4 * defaultPar.rho;
    case "NLM"
        defaultPar.xi = 1e-7;
        defaultPar.rho = 0.002;
        defaultPar.lambda = 0.01 * defaultPar.rho;
end
addOptional(p,'rho',defaultPar.rho,@isnumeric);
addOptional(p,'lambda',defaultPar.lambda,@isnumeric);
addOptional(p,'xi',defaultPar.xi,@isnumeric);
addOptional(p,'init_rho',0,@isnumeric);
addOptional(p,'init_xi',0,@isnumeric);
addOptional(p,'ShowInfo',false,@islogical);

parse(p,varargin{:});
par = p.Results;


rho = par.rho;
lambda = par.lambda;
xi = par.xi;
init_rho = par.init_rho * ones(R,1);
init_xi = par.init_xi;
X_true = par.X_true;
S_true = par.S_true;
C_true = par.C_true;
BuildingMask = par.BuildingMask;
max_iternum = par.max_iter;
KEEP_W = ~par.NLM_MEMORY_EFF;

th = par.threshold;
ShowInfo = par.ShowInfo; 


%% Init
start_pnp = tic;

[I,J,K] = size(X_observed);
x_max = max(X_observed(:));
X_observed = X_observed ./ x_max ./ 10;

if all(BuildingMask==0) % all zero, no building info
    [S_est, C_est, R_est] = initial_pnp(X_observed, Omega, R, init_rho, init_xi);
    R = R_est;
    realOmega = Omega;
    ind_realO = find(realOmega);
else
    Omega_aug = zeros(I,J,2);
    Omega_aug(:,:,1) = Omega; Omega_aug(:,:,2) = BuildingMask;
    [S_est, C_est, R_est] = initial_pnp(X_observed, Omega_aug, R, init_rho, init_xi);
    R = R_est;
    ind_O = find(Omega);
    ind_B = find(BuildingMask);
    ind_realO = setdiff(ind_O,ind_B);
    realOmega = zeros(I,J);
    realOmega(ind_realO) = 1;
end

S_mat = reshape(S_est,[I*J,R]);
X_mat = S_mat*C_est';
X_est = reshape(X_mat,[I,J,K]);

W = zeros(I,J,R);
X_last = 0;
S_last = 0;
Z_last = 0;
W_last = 0;
dist = zeros(max_iternum,1);
dist_all = zeros(max_iternum,1);

rhoList = rho * ones(R,1);

%% LaPnP
for i = 1:max_iternum
    switch Denoiser
        case "BM3D"
            [Z_est,~] = BM3D_denoising(S_est,W,rho,lambda); % I * J * R
        case "DPIRNet"
            Z_est = dpir_denoising(S_est,W,rho,lambda,BuildingMask);
        case "NLM"
            if KEEP_W
                if i <= 10
                    [Z_est,denoiser_W] = NLM_denoising(S_est,W,rho,lambda);
                else
                   Z_est = NLM_denoising(S_est,W,rho,lambda,denoiser_W);
                end
            else
                Z_est = NLM_denoising(S_est,W,rho,lambda);
            end
    end

    if ShowInfo
        figure(4)
        contourSC(0,S_est(:,:,1),Z_est(:,:,1),S_est(:,:,2),Z_est(:,:,2),S_est(:,:,3),Z_est(:,:,3));
    end

    [X_est, S_est, C_est] = prox_f(X_observed, Omega, Z_est - W, ...
            S_est, C_est, rhoList, xi);


    W = W + S_est - Z_est;
    
    dist(i) = norm(X_est(:) - X_last(:));
    dist_all(i) = norm(S_est(:) - S_last(:)) + norm(Z_est(:)-Z_last(:)) + norm(W(:)-W_last(:));
    dist(i) = dist(i)/norm(X_last(:));

    if i >= 2 && dist_all(i) >= 0.95 * dist_all(i-1) ...
             && ~strcmp(par.Denoiser,"NLM")
        rho = rho*1.1;
    end

    if dist(i) < th && i >=5
        break
    end

    X_last = X_est;
    S_last = S_est;
    Z_last = Z_est;
    W_last = W;
  
end

X_est = X_est .* 10 .* x_max;
S_est = S_est .* 10 .* x_max;
time_est = toc(start_pnp);
end