clear all;

% functions for admm-pnp
addpath(('f_data'))
addpath(genpath('f_alg'))
addpath('f_evaluation')

% Load BM3D denoiser
if ismac
    addpath(genpath(fullfile(pwd, 'Denoisers', 'bm3d_matlab_package_4.0.3')));
else
    addpath(genpath(fullfile(pwd, 'Denoisers', 'BM3D')));
end


% Load DRUnet, utils has been changed to utils1
addpath(fullfile(pwd, 'Denoisers', 'DPIR-master'));
py_path = { ...
    fullfile(pwd, 'Denoisers', 'DPIR-master'), ...
    fullfile(pwd, 'Denoisers', 'DPIR-master', 'utils1'), ...
    fullfile(pwd, 'Denoisers', 'DPIR-master', 'model_zoo') ...
};

for i = 1:length(py_path)
    if count(py.sys.path,py_path{i}) == 0
        insert(py.sys.path,int32(0),py_path{i});
    end
end
dpirModel = py.importlib.import_module('dpir_denoising_interface');
py.importlib.reload(dpirModel);

%% user define
map_size = [51,51,32];
R = 3;
shadow_sigma = 6;
d_corr =  50;
use_snr = false;
snr = 20;
sampleRate = 0.1;

%% data generation
[X_true,Sc,C] = smap_gen(map_size, R, shadow_sigma, d_corr, 'num_peaks_per_psd',3,'strictly_separable',false);

[I,J,K] = size(X_true);

% Add noise
if use_snr
    Ps = norm(X_true(:))^2;
    Pn = Ps*10^(-snr/10);

    sn = sqrt(Pn/numel(X_true));
    Noise = (sn*normrnd(0,1,[I,J,K]));
    X_noisy = max(X_true + Noise,0);
else
    X_noisy = X_true;
end

% Samples of X, random fibers
Omega = zeros(map_size([1,2]));
if ~exist('BuildingMask','var')
    BuildingMask = zeros(map_size([1,2]));
    mask_index = [];
else
    mask_index = find(BuildingMask);
end

noMask_index = find(1-BuildingMask);
num_samples = ceil(sampleRate*length(noMask_index));
ind_O = randperm(   length(noMask_index), num_samples   )';
ind_realO = noMask_index(ind_O);
Omega(ind_realO) = 1;

Omega_aug = Omega;
Omega_aug(:,:,2) = BuildingMask;

X_observed = X_noisy;
for k = 1:map_size(3)
    X_observed(:,:,k) = X_observed(:,:,k).*Omega.*(1-BuildingMask);
end

%% recover X

% -------- Parameters set for synthetic data without separable conditions, otherwise use the builtin default setting ----------
rho = 0.1 * max(1,(260.1/length(ind_realO)));
xi = 1e-3 * R/6;

% -------- bm3d ----------
lambda = 10 *rho;
[X_bm3d, time_bm3d, S_bm3d, C_bm3d] = sc_pnp(X_observed, Omega, R, "BM3D", ...
                                            'BuildingMask',BuildingMask,...
                                            'rho',rho,...
                                            'lambda',lambda,....
                                            'xi',xi);
X_bm3d = X_bm3d.*(1-BuildingMask);


% -------- dpir ----------
lambda = 2.5e-4 *rho;
[X_dpir, time_dpir, S_dpir, C_dpir] = sc_pnp(X_observed, Omega, R, "DPIRNet",...
                                            'BuildingMask',BuildingMask,...
                                            'rho',rho,...
                                            'lambda',lambda,....
                                            'xi',xi);
X_dpir = X_dpir.*(1-BuildingMask);

% -------- nlm ----------
lambda = 0.01 *rho;
[X_nlm, time_nlm, S_nlm, C_nlm] = sc_pnp(X_observed, Omega, R, "NLM",...
                                        'BuildingMask',BuildingMask,...
                                        'rho',rho,...
                                        'xi',xi,...
                                        'lambda',lambda);
X_nlm = X_nlm.*(1-BuildingMask);


%% evaluation
kk = 10; figure;
contourSC(BuildingMask,(1-BuildingMask).*X_true(:,:,kk),"Ground truth",...
          (1-BuildingMask).*X_observed(:,:,kk),"Observed",...
          (1-BuildingMask).*X_nlm(:,:,kk),"LaPNP-NLM",...
          (1-BuildingMask).*X_dpir(:,:,kk),"LaPnP-DRUnet",...
          (1-BuildingMask).*X_bm3d(:,:,kk),"LaPnP-BM3D")

fprintf('rse between X_observed and X_true: %.3f\n', norm(X_observed(:)-X_true(:))^2/norm(X_true(:))^2);
fprintf('rse between X_nlm and X_true: %.3f\n', norm(X_nlm(:)-X_true(:))^2/norm(X_true(:))^2);
fprintf('rse between X_dpir and X_true: %.3f\n', norm(X_dpir(:)-X_true(:))^2/norm(X_true(:))^2);
fprintf('rse between X_bm3d and X_true: %.3f\n', norm(X_bm3d(:)-X_true(:))^2/norm(X_true(:))^2);

fprintf('ssim between X_observed and X_true: %.4f\n', ssim_score_4sc(X_observed,X_true) );
fprintf('ssim between X_nlm and X_true: %.4f\n', ssim_score_4sc(X_nlm,X_true) );
fprintf('ssim between X_dpir and X_true: %.4f\n', ssim_score_4sc(X_dpir,X_true));
fprintf('ssim between X_bm3d and X_true: %.4f\n', ssim_score_4sc(X_bm3d,X_true) );