clear all;

% functions for admm-pnp
addpath(('f_data'))
addpath(genpath('f_alg'))
addpath('f_evaluation')

% Load BM3D denoiser
addpath(genpath(fullfile(pwd, 'Denoisers', 'BM3D')));

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
use_snr = false;
snr = 20;
sampleRate = 0.05;

%% data generation
[X_true,BuildingMask] = smap_gen_Mannheim();

[I,J,K] = size(X_true);
map_size = [I,J,K];
R = 7;

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
Omega = zeros(I,J);
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
% -------- bm3d ----------
[X_bm3d, time_bm3d, S_bm3d, C_bm3d] = sc_pnp(X_observed, Omega, R, "BM3D", ...
                                            'BuildingMask',BuildingMask);
X_bm3d = X_bm3d.*(1-BuildingMask);

% -------- dpir ----------
[X_dpir, time_dpir, S_dpir, C_dpir] = sc_pnp(X_observed, Omega, R, "DPIRNet",...
                                            'BuildingMask',BuildingMask);
X_dpir = X_dpir.*(1-BuildingMask);

% -------- nlm ----------
[X_nlm, time_nlm, S_nlm, C_nlm] = sc_pnp(X_observed, Omega, R, "NLM",...
                                        'BuildingMask',BuildingMask);
X_nlm = X_nlm.*(1-BuildingMask);


%% evaluation
contourSC(kron(ones(K,1),BuildingMask),reshape(permute(X_true,[1,3,2]),[I*K,J]),"Ground truth",...
          reshape(permute(X_observed,[1,3,2]),[I*K,J]),"Observed",...
          reshape(permute(X_nlm,[1,3,2]),[I*K,J]),"LaPNP-NLM",...
          reshape(permute(X_dpir,[1,3,2]),[I*K,J]),"LaPnP-DRUnet",...
          reshape(permute(X_bm3d,[1,3,2]),[I*K,J]),"LaPnP-BM3D")


fprintf('rse between X_observed and X_true: %.3f\n', norm(X_observed(:)-X_true(:))^2/norm(X_true(:))^2);
fprintf('rse between X_nlm and X_true: %.3f\n', norm(X_nlm(:)-X_true(:))^2/norm(X_true(:))^2);
fprintf('rse between X_dpir and X_true: %.3f\n', norm(X_dpir(:)-X_true(:))^2/norm(X_true(:))^2);
fprintf('rse between X_bm3d and X_true: %.3f\n', norm(X_bm3d(:)-X_true(:))^2/norm(X_true(:))^2);

fprintf('ssim between X_observed and X_true: %.4f\n', ssim_score_4sc(X_observed,X_true) );
fprintf('ssim between X_nlm and X_true: %.4f\n', ssim_score_4sc(X_nlm,X_true) );
fprintf('ssim between X_dpir and X_true: %.4f\n', ssim_score_4sc(X_dpir,X_true));
fprintf('ssim between X_bm3d and X_true: %.4f\n', ssim_score_4sc(X_bm3d,X_true) );