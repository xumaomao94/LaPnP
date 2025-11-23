function [y_est, matches] = BM4D(z, sigma_psd, profile, stage_arg, blockmatches)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  BM4D is an algorithm for attenuation of additive spatially correlated
%  stationary (aka colored) Gaussian noise in volumetric data.
%
%
%  FUNCTION INTERFACE:
%
%  y_est = BM4D(z, sigma_psd, profile)
%
%  INPUT ARGUMENTS:
%
%  -- required --
%
%         'z' : noisy image (3-D array)
%  'sigma_psd' : noise power spectral density (3-D nonnegative array)
%               OR
%               noise STD
%
% -- optional --
%
%   'profile' : 'np' --> Normal Profile (default)
%               'refilter' --> Apply refiltering
%               OR
%               a BM4DProfile object specifying the parameters
%               some other premade profiles also included from the previous versions
%               in BM4DProfile.m
%
%   'stage_arg' : Determines whether to perform hard-thresholding or wiener filtering.
%                 either BM4DProfile.HARD_THRESHOLDING, BM4DProfile.ALL_STAGES or an estimate
%                  of the noise-free image.
%                    - BM4DProfile.ALL_STAGES: Perform both.
%                    - BM4DProfile.HARD_THRESHOLDING: Perform hard-thresholding only.
%                    - ndarray, size of z: Perform Wiener Filtering with stage_arg as pilot.
%   'blockmatches' : Blockmatch data pair (cell), one element for HT, one
%                   for Wie. 
%                   - false: don't store block matches
%                   - true: store block matches, return value will be
%                   [y_est, matches]
%                   - blockmatches struct returned by previous application,
%                   will be used instead of matching on the input
%
%  OUTPUT:
%      'y_est'  denoised image  (equal size to z)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2006-2021 Tampere University.
% All rights reserved.
% This work (software, material, and documentation) shall only
% be used for nonprofit noncommercial purposes.
% Any unauthorized use of this work for commercial or for-profit purposes
% is prohibited.
%
% AUTHORS:
%     Y. Mäkinen, L. Azzari, K. Dabov, A. Foi
%     email: ymir.makinen@tuni.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('profile','var')
    profile         = 'np'; %% default profile
end

if isa(profile, 'string') || isa(profile, 'char')
    profile = BM4DProfile(profile);
end

if ~exist('blockmatches','var')
    blockmatches = {false, false};
end
if numel(blockmatches) ~= 2
    disp('Error: two values are needed for block-match data! (use "false" to ignore)');
end

if ~exist('stage_arg','var')
    stage_arg = profile.ALL_STAGES;  % By default, do both HT and Wie
elseif stage_arg == profile.WIENER_FILTERING
    disp('Error: If you wish to only do wiener filtering, pass the estimate y_hat instead of the WIENER_FILTERING value!')
    return
elseif isa(stage_arg, 'float') || isa(stage_arg, 'double')
    
    if numel(size(stage_arg)) < 2
        disp('Error: stage_arg must be either stage value from BM4DProfile or an estimate y_hat!')
        return
    end
    
    % Presume that stage_arg is an estimate for wiener.
    y_hat = stage_arg;
    stage_arg = profile.WIENER_FILTERING;

elseif ~isa(stage_arg, 'int8')
    disp('Error: stage_arg must be either stage value from BM4DProfile or an estimate y_hat!')
    return
end

for i=1:3
    if size(z, i) < profile.N1(i) || size(z, i) < profile.N1_wiener(i)
        disp('Error: Image cannot be smaller than block size!')

    if all(profile.N1(1:2) <= size(z, 1:2)) && all(profile.N1_wiener(1:2) <= size(z, 1:2))
       disp('(you may be using a 3-D filtering profile for 2-D data)'); 
    end

    return
    end
end
% Define maximum pad size: pad size should be at least
% half the size of the correlation kernel, but needn't be larger
% (although it can be)
% To be sure, we pad total of the image size, but if the
% kernel size is approximately known, some computation time
% may be saved by specifying it in the profile.
if profile.max_pad_size(1) == -1
    pad_size = ceil(size(z)/2);
else
    pad_size = profile.max_pad_size;
end

% Conventional mode
if sum(profile.Nf) == 0
    profile.Nf = [min(size(z, 1), 16), min(size(z, 2), 16), min(size(z, 3), 16)];
    profile.Kin = 0;
    profile.gamma = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

single_dim_psd = false;
% Check if a sigma was passed instead of a PSD
if(numel(sigma_psd) == 1)
    single_dim_psd = true;
    sigma_psd = ones(size(z)) * sigma_psd^2 * numel(z);
end

% Get relevant info from the sigma_psd, including lambda and mu.
[sigma_psd2, psd_blur, psd_k, profile] = process_psd(sigma_psd, z, ...
    single_dim_psd, pad_size, profile);

bms_ht_out = false;
bms_wiener_out = false;

collect_bm_ht = int32(-1);
collect_bm_wiener = int32(-1);
bms_in_ht = false;
bms_in_wiener = false;

if isa(blockmatches{1}, 'struct')
    collect_bm_ht = int32(1);
    bms_in_ht = blockmatches{1};
elseif blockmatches{1}
    collect_bm_ht = int32(0);
end

if isa(blockmatches{2}, 'struct')
    collect_bm_wiener = int32(1);
    bms_in_wiener = blockmatches{2};
elseif blockmatches{2}
    collect_bm_wiener = int32(0);
end

is_real = isreal(z);

mu2DC = profile.mu2;
mu2_reDC = profile.mu2_re;

profile.splitComp(1) = false;

size_z = ones(1, 3);
size_z(1:numel(size(z))) = size(z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 1. Produce the basic estimate by HT filtering
%%%%

do_sharpen = profile.sharpen_alpha ~= 1 || profile.sharpen_alpha_3d ~= 1;
do_sharpen_aggregation = do_sharpen;

if bitand(stage_arg, BM4DProfile.HARD_THRESHOLDING)
    
    qshifts = get_shift_params(profile.N1, profile.Nstep);

    [t_forward, t_inverse, hadper_trans_single_den, ...
        inverse_hadper_trans_single_den, Wwin3D] = get_transforms(profile, true);

    [y_hat, bms_ht_out] = bm4dmx(single(z), hadper_trans_single_den, ...
        inverse_hadper_trans_single_den, int32(profile.Nstep), int32(profile.N1), ... 
        int32(profile.N2), ...
        single(profile.lambda_thr), single(profile.tau_match*(max(z(:)) - min(z(:)))^2), ...
        int32(profile.Ns), t_forward, t_inverse,...
        single(Wwin3D), single(psd_blur), ...
        int32(profile.Nf), single(profile.gamma), int32(profile.Kin), ...
        int8(profile.splitComp), qshifts, collect_bm_ht, bms_in_ht, int32(0), [], ...
        do_sharpen, do_sharpen_aggregation, single(1/profile.sharpen_alpha), single(1/profile.sharpen_alpha_3d), profile.num_threads);


    % Re-filter    
    if (profile.denoise_residual)

       [remains, remains_PSD] = get_filtered_residual(z, y_hat, sigma_psd2, pad_size, ...
           profile.residual_thr, single_dim_psd);
       remains_PSD = process_psd_for_nf(remains_PSD, psd_k, profile);
        
        % Skip refiltering if there is a zero sigma_psd
        if(min(max(max(remains_PSD, [], 1), [], 2)) > 1e-5)
            
            z_2 = single(y_hat + remains);
            
            % Re-filter
            tau_scale = max(z_2(:)) - min(z_2(:));
            [y_hat, bms_ht_out] = bm4dmx(z_2,hadper_trans_single_den, ...
                    inverse_hadper_trans_single_den, int32(profile.Nstep), int32(profile.N1), ... 
                    int32(profile.N2), ...
                    single(profile.lambda_thr_re), single(profile.tau_match*tau_scale^2), ...
                    int32(profile.Ns), t_forward, t_inverse,...
                    single(Wwin3D), single(remains_PSD), ...
                    int32(profile.Nf), single(profile.gamma), int32(profile.Kin), ...
                    int8(profile.splitComp), qshifts, collect_bm_ht, bms_in_ht, int32(0), [], ...
                    do_sharpen, do_sharpen_aggregation, single(1/profile.sharpen_alpha), single(1/profile.sharpen_alpha_3d), profile.num_threads);

        end


    end

    if profile.print_info
        disp('Hard-thresholding phase completed')
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 2. Produce the final estimate by Wiener filtering (using the
%%%%  hard-thresholding initial estimate)
%%%

if bitand(stage_arg, BM4DProfile.WIENER_FILTERING)

    qshifts = get_shift_params(profile.N1_wiener, profile.Nstep_wiener);

    [t_forward, t_inverse, hadper_trans_single_den, ...
        inverse_hadper_trans_single_den, Wwin3D] = get_transforms(profile, false);
    % Wiener filtering
    
    tau_scale = max(y_hat(:)) - min(y_hat(:));
    [y_est, bms_wiener_out] = bm4dmx(single(z),hadper_trans_single_den, ...
                    inverse_hadper_trans_single_den, int32(profile.Nstep_wiener), ...
                    int32(profile.N1_wiener), ... 
                    int32(profile.N2_wiener), ...
                    single([profile.mu2, mu2DC]), single(profile.tau_match_wiener*tau_scale^2), ...
                    int32(profile.Ns_wiener), t_forward, t_inverse,...
                    single(Wwin3D), single(psd_blur), ...
                    int32(profile.Nf), single(0), int32(profile.Kin), ...
                    int8(profile.splitComp), qshifts, collect_bm_wiener, bms_in_wiener, ...
                    int32(1), single(y_hat), false, false, 1, 1, profile.num_threads);



    if (profile.denoise_residual)

        [remains, remains_PSD] = get_filtered_residual(z, y_est, sigma_psd2, pad_size, ...
            profile.residual_thr, single_dim_psd);
        remains_PSD = process_psd_for_nf(remains_PSD, psd_k, profile);
        
        if(min(max(max(remains_PSD, [], 1), [], 2)) > 1e-5) 

            % Re-filter
            
            tau_scale = max(y_est(:)) - min(y_est(:));
            [y_est, bms_wiener_out] = bm4dmx(single(y_est + remains),hadper_trans_single_den, ...
                    inverse_hadper_trans_single_den, int32(profile.Nstep_wiener), ...
                    int32(profile.N1_wiener), ... 
                    int32(profile.N2_wiener), ...
                    single([profile.mu2_re, mu2_reDC]), single(profile.tau_match_wiener*tau_scale^2), ...
                    int32(profile.Ns_wiener), t_forward, t_inverse,...
                    single(Wwin3D), single(remains_PSD), ...
                    int32(profile.Nf), single(0), int32(profile.Kin), ...
                    int8(profile.splitComp), qshifts, collect_bm_wiener, bms_in_wiener, ...
                    int32(1), single(y_est), false, false, 1, 1, profile.num_threads);
        end

    end    

    if profile.print_info
        disp('Wiener phase completed')
    end
else
    y_est = y_hat;
end

y_est = double(y_est);
if collect_bm_ht == 0 || collect_bm_wiener == 0
    matches = {false, false};
    if collect_bm_ht == 0
        matches{1} = bms_ht_out;
    end
    if collect_bm_wiener == 0
        matches{2} = bms_wiener_out;
    end
end

if is_real
    % Input was real, output should be real
    y_est = real(y_est);
end

return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Get residual, filtered by global FFT HT
function [remains, remains_PSD] = get_filtered_residual(z, y_hat, sigma_psd, pad_size, ...
    residual_thr, single_dim_psd)
    
    resid = fftn(padarray(z - double(y_hat), pad_size, 'both'));

    % Convolve & max filter the thresholded map for extra coefficients
    ksz = ceil(size(resid) / 150);
    ksz = ksz + 1 - mod(ksz, 2);
    k=fspecial('gaussian',[ksz(1), 1], size(resid, 1) / 500) * ...
        fspecial('gaussian',[1, ksz(2)], size(resid, 2) / 500);
    if numel(ksz) > 2
        k=repmat(k, [1, 1, ksz(3)]) .* ...
            permute(fspecial('gaussian',[ksz(3), 1], size(resid, 1) / 500), [2, 3, 1]);
    end
    if numel(pad_size) == 2
        pad_size(3) = 0;
    end
    P = sigma_psd;

    cc=convn(padarray(abs(resid) > residual_thr .* sqrt(P), (size(k) - 1)/2, 'circular'), k, 'valid');

    % Threshold mask
    msk = (cc > 0.01);
    
    % Residual + sigma_psd
    remains = ifftn((resid) .* msk);
    remains_PSD = P .* msk;

    % Crop the pad off
    remains = remains(1 + pad_size(1) : end - pad_size(1),...
                      1 + pad_size(2) : end - pad_size(2), ...
                      1 + pad_size(3) : end - pad_size(3));
   
    % Also from the sigma_psd             
    temp_kernel = real(fftshift(ifftn(sqrt(remains_PSD / (numel(z))))));
    temp_kernel = temp_kernel(1+pad_size(1):end-pad_size(1), 1+pad_size(2):end-pad_size(2),  1+pad_size(3):end-pad_size(3));
    remains_PSD = abs(fftn(temp_kernel, size(z))).^2 * numel(z);

end

function [Tforward, Tinverse] = get_transf_matrix(N, transform_type, dec_levels)
%
% Create forward and inverse transform matrices, which allow for perfect
% reconstruction. The forward transform matrix is normalized so that the
% l2-norm of each basis element is 1.
%
% [Tforward, Tinverse] = get_transf_matrix (N, transform_type, dec_levels)
%
%  INPUTS:
%
%   N               --> Size of the transform (for wavelets, must be 2^K)
%
%   transform_type  --> 'dct', 'dst', 'hadamard', or anything that is
%                       listed by 'help wfilters' (bi-orthogonal wavelets)
%                       'DCrand' -- an orthonormal transform with a DC and all
%                       the other basis elements of random nature
%
%   dec_levels      --> If a wavelet transform is generated, this is the
%                       desired decomposition level. Must be in the
%                       range [0, log2(N)-1], where "0" implies
%                       full decomposition.
%
%  OUTPUTS:
%
%   Tforward        --> (N x N) Forward transform matrix
%
%   Tinverse        --> (N x N) Inverse transform matrix
%

if ~exist('dec_levels','var')
    dec_levels = 0;
end

if N == 1
    Tforward = 1;
    Tinverse = 1;
    return;
elseif strcmp(transform_type, 'hadamard') == 1
    Tforward = hadamard(N)/sqrt(N);
    Tinverse = Tforward;
    return;
elseif (N == 8) && strcmp(transform_type, 'bior1.5')==1 % hardcoded transform so that the wavelet toolbox is not needed to generate it
    Tforward =[ 0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110
               -0.225454819240296  -0.461645582253923  -0.461645582253923  -0.225454819240296   0.225454819240296   0.461645582253923   0.461645582253923   0.225454819240296
                0.569359398342840   0.402347308162280  -0.402347308162280  -0.569359398342840  -0.083506045090280   0.083506045090280  -0.083506045090280   0.083506045090280
               -0.083506045090280   0.083506045090280  -0.083506045090280   0.083506045090280   0.569359398342840   0.402347308162280  -0.402347308162280  -0.569359398342840
                0.707106781186550  -0.707106781186550                   0                   0                   0                   0                   0                   0
                                0                   0   0.707106781186550  -0.707106781186550                   0                   0                   0                   0
                                0                   0                   0                   0   0.707106781186550  -0.707106781186550                   0                   0
                                0                   0                   0                   0                   0                   0   0.707106781186550  -0.707106781186550];
elseif strcmp(transform_type, 'dct') == 1
    Tinverse = cos(pi*(0:N-1).*(2*(0:N-1)'+1)/(2*N))./sqrt(N./[1 2*ones(1,N-1)]);
    Tforward = Tinverse';
    return;
elseif strcmp(transform_type, 'dst') == 1
    Tforward = sin(pi*(1:N).*(1:N)'/(N+1))*sqrt(2/(N+1));
    Tinverse = Tforward';
    return;
elseif strcmp(transform_type, 'fft') == 1
    Tforward    = fft(eye(N)) / sqrt(N);
    Tinverse    = Tforward';
    return;
elseif strcmp(transform_type, 'DCrand') == 1
    x = randn(N); x(1:end,1) = 1; [Q,~] = qr(x);
    if (Q(1) < 0)
        Q = -Q;
    end
    Tforward = Q';
else %% a wavelet decomposition supported by 'wavedec'
    %%% Set periodic boundary conditions, to preserve bi-orthogonality
    dwtmode('per','nodisp');

    Tforward = zeros(N,N);
    for i = 1:N
        Tforward(:,i)=wavedec(circshift([1 zeros(1,N-1)],[dec_levels i-1]), log2(N), transform_type);  %% construct transform matrix
    end
end

%%% Normalize the basis elements
if ~((N == 8) && strcmp(transform_type, 'bior1.5')==1)
    Tforward = (Tforward' * diag(sqrt(1./sum(Tforward.^2,2))))';
end

%%% Compute the inverse transform matrix
Tinverse = inv(Tforward);

end


function [psd] = process_psd_for_nf(sigma_psd, psd_k, profile)
    if profile.Nf == 0
        psd = sigma_psd;
        return;
    end
    
    % Reduce PSD size to start with
    max_ratio = 16;
    sigma_psd_copy = sigma_psd;
    single_kernel_x = ones(3, 1, 1) / 3;
    single_kernel_y = ones(1, 3, 1) / 3;
    single_kernel_z = ones(1, 1, 3) / 3;
    single_kernels = {single_kernel_x, single_kernel_y, single_kernel_z};
    for k = 1:3
        orig_ratio = size(sigma_psd, k) / profile.Nf(k);
        ratio = orig_ratio;
        single_kernel = single_kernels{k};
        while ratio > max_ratio
            mid_corr = convn(padarray(sigma_psd_copy, double([k==1, k==2, k==3]), 'circular'), ...
                             single_kernel, 'valid');
            if k == 1             
                sigma_psd_copy = mid_corr(2:3:end, :, :);
            elseif k == 2
                sigma_psd_copy = mid_corr(:, 2:3:end, :);
            else
                sigma_psd_copy = mid_corr(:, 2:3:end, 2:3:end);
            end
            ratio = size(sigma_psd_copy, k) / profile.Nf(k);
        end
    
        % Scale PSD because the binary expects it to be scaled by size
        sigma_psd_copy = sigma_psd_copy .* (ratio / orig_ratio);
        
    end
    if ~isempty(psd_k)
        sigma_psd_copy = convn(padarray(sigma_psd_copy, ...
                                        (size(psd_k) - 1) / 2, 'circular'), ...
                               psd_k, 'valid');
    end
    
    psd = sigma_psd_copy;

end

function [lambda, wielambdasq, lambda2, wielambdasq2] = estimate_psd_parameters(PSD65_full)
% Estimate parameters based on the sigma_psd

% Get the optimal parameters and matching features for a bunch of PSDs
load('param_matching_data.mat', 'features', 'maxes');

sz = 65;
indices_to_take = [1:2:10 12:5:32];

lambda = [];
wielambdasq = [];
lambda2 = [];
wielambdasq2 = [];

% Get separate parameters for each sigma_psd provided
for psd_num = 1:size(PSD65_full, 3)
    PSD65 = PSD65_full(:, :, psd_num);
    PSD65 = fftshift(PSD65);
    
    % Get features for this sigma_psd
    pcaxa = get_features(PSD65, sz, indices_to_take);

    % Calculate distances to other PSDs
    mm = mean(features, 2);
    centered_features = features - mm;
    corr_matx = centered_features * centered_features';
    corr_matx = corr_matx / 500;
    
    centered_pcax = pcaxa' - mm;
    
    [u, s, ~] = svd(corr_matx);
    centered_features = u * centered_features;
    centered_pcax = u * centered_pcax;
    centered_features = centered_features .* sqrt(diag(s));
    centered_pcax = centered_pcax .* sqrt(diag(s));

    diff_pcax = sqrt(sum(abs(centered_features - centered_pcax).^2, 1));

    % Take only smallest->best x %
    [~, dff_I] = sort(diff_pcax);

    %  Take 20 most similar PSDs into consideration
    count = 20;
    diff_indices = dff_I(1:count);

    % Invert, smaller -> bigger weight
    diff_inv = 1 ./ (diff_pcax + eps);
    diff_inv = diff_inv(diff_indices) ./ sum(diff_inv(diff_indices));

    % Weight
    param_idxs = sum(diff_inv .* maxes(diff_indices, :)', 2);

    lambdas = 2.5:0.1:5;
    wielambdasqs = 0.2:0.2:6;
    
    % Get parameters from indices - 
    % Interpolate lambdas and mu^2s from the list
    for ix = [1, 3]
        param_idx = max(1, param_idxs(ix));
        param_idx2 = max(1, param_idxs(ix+1));

        l1 = lambdas(floor(param_idx));
        l2 = lambdas(min(ceil(param_idx), numel(lambdas)));

        w1 = wielambdasqs(floor(param_idx2));
        w2 = wielambdasqs(min(ceil(param_idx2), numel(wielambdasqs)));

        param_smooth = param_idx - floor(param_idx);
        param_smooth2 = param_idx2 - floor(param_idx2);

        if ix == 1
            lambda = [lambda, l2 * param_smooth + l1 * (1 - param_smooth)];
            wielambdasq = [wielambdasq, w2 * param_smooth2 + w1 * (1 - param_smooth2)];
        elseif ix == 3
            lambda2 = [lambda2, l2 * param_smooth + l1 * (1 - param_smooth)];
            wielambdasq2 = [wielambdasq2, w2 * param_smooth2 + w1 * (1 - param_smooth2)];
        end

    end
end
end

% Calculate features for a sigma_psd from integrals
function f = get_features(sigma_psd, sz, indices_to_take)

    [I_rot, I_rot2] = pcax(sigma_psd);
    f1 = zeros(1, numel(indices_to_take));
    f2 = f1;
    
    % Extract features for a sigma_psd
    for ii = 1:numel(indices_to_take)
        rang = indices_to_take(ii);
        if ii > 1
            rang = indices_to_take(ii-1)+1:rang;
        end
        f1(ii) = sum(I_rot(ceil(sz/2) + rang - 1)) / numel(rang);
        f2(ii) = sum(I_rot2(ceil(sz/2) + rang - 1)) / numel(rang);
    end
    
    f = [f1 f2];
end

% Calculate integrals along principal axes of the sigma_psd
function [I_rot, I_rot2] = pcax(sigma_psd)

N=size(sigma_psd,1);
[G2, G1]=meshgrid(1:N,1:N);

trapz2D=@(G2,G1,sigma_psd) trapz(G1(:,1),trapz2(G2,sigma_psd,2),1);

Pn=sigma_psd/trapz2D(G2,G1,sigma_psd);

m2=trapz2D(G2,G1,Pn.*G2);
m1=trapz2D(G2,G1,Pn.*G1);
C=zeros(2);
O1=[2 1 1 0];
O2=[0 1 1 2];
for jj=[1 2 4]
    C(jj)=squeeze(trapz2D(G2,G1,Pn.*(G2-m2).^O1(jj).*(G1-m1).^O2(jj)));
end
C(3)=C(2);


[U, ~, ~]=svd(C);

N3 = 3 * N;
[G13N, G23N]=ndgrid((1:N3)-(N3+1)/2, (1:N3)-(N3+1)/2);

% Rotate PSDs and calculate integrals along the rotated PSDs
theta = angle(U(1, 1) + 1i * U(1, 2));
G2rot=G23N(N+1:2*N,N+1:2*N)*cos(theta)-G13N(N+1:2*N,N+1:2*N)*sin(theta);
G1rot=G13N(N+1:2*N,N+1:2*N)*cos(theta)+G23N(N+1:2*N,N+1:2*N)*sin(theta);
P_rot_handle = griddedInterpolant(G13N,G23N,repmat(sigma_psd,[3,3]),'linear','nearest');
P_rot = P_rot_handle(G1rot,G2rot);
I_rot = trapz2(G1, P_rot, 1);

theta2 = angle(U(2, 1) + 1i * U(2, 2));
G2rot=G23N(N+1:2*N,N+1:2*N)*cos(theta2)-G13N(N+1:2*N,N+1:2*N)*sin(theta2);
G1rot=G13N(N+1:2*N,N+1:2*N)*cos(theta2)+G23N(N+1:2*N,N+1:2*N)*sin(theta2);
P_rot2 = P_rot_handle(G1rot,G2rot);
I_rot2 = trapz2(G1, P_rot2, 1);
end

function I = trapz2(X,Y,dimm)

if dimm==2
I=sum((Y(:,2:end)+Y(:,1:end-1))/2.*(X(:,2:end)-X(:,1:end-1)),2);
else
I=sum((Y(2:end,:)+Y(1:end-1,:))/2.*(X(2:end,:)-X(1:end-1,:)),1);
end

end

% Process sigma_psd, get parameters if needed
function [sigma_psd2, psd_blur, psd_k, profile] = process_psd(sigma_psd, z, single_dim_psd, ...
    pad_size, profile)

% Get auto params if there is a relevant parameter to be calculated.
auto_params = profile.lambda_thr == profile.NO_VALUE || ...
              profile.mu2 == profile.NO_VALUE || ...
              (profile.denoise_residual && ...
                (profile.lambda_thr_re == profile.NO_VALUE || ...
                profile.mu2_re == profile.NO_VALUE) ...
              );


% Calculate the correlation kernel from the sigma_psd in order
% to resize the sigma_psd. (skip if we are not resizing the sigma_psd)
if (profile.denoise_residual || auto_params)
    temp_kernel = real(fftshift(ifftn(sqrt(sigma_psd / numel(sigma_psd)))));
end

% We need a bigger sigma_psd if we are denoising residuals
if profile.denoise_residual && pad_size(1)
    extended_size = size(z) + pad_size*2;
    % bigger sigma_psd
    sigma_psd2 = abs(fftn(temp_kernel, extended_size)).^2 * numel(z);
else
    sigma_psd2 = sigma_psd;
end

if auto_params && ~single_dim_psd
    % Estimate parameters based on the sigma_psd
    
    minus_size = ceil((size(z) - 65) / 2);
    temp_kernel65 = temp_kernel;
    if minus_size(1) > 0
        temp_kernel65 = temp_kernel65(1+minus_size(1):minus_size(1)+65, :, :);
    end
    if minus_size(2) > 0
        temp_kernel65 = temp_kernel65(:, 1+minus_size(2):minus_size(2)+65, :);
    end
    if numel(minus_size) > 3 && minus_size(3) > 0
        temp_kernel65 = temp_kernel65(:, :, 1+minus_size(3):minus_size(3)+65);
    end
    temp_kernel65 = temp_kernel65 / norm(temp_kernel65(:));
    if size(temp_kernel65, 3) > 1
        PSD65 = abs(fftn(temp_kernel65, [65, 65, 65])).^2 * 65 * 65;
        % Reduce PSD to two dimensions to use BM3D parameter estimation
        PSD65 = project_psd_for_param_est(fftshift(PSD65));
    else
        PSD65 = abs(fftn(temp_kernel65, [65, 65])).^2 * 65 * 65;
    end
    PSD65 = PSD65 / mean(PSD65(:)) * 65 * 65;

    % Parameter estimation
    [lambda_thr, mu2, lambda_thr_re, mu2_re] = estimate_psd_parameters(PSD65);
else
    % For white noise, the result of the estimation is this.
    % No need to create a sigma_psd and run the script just for that.
    lambda_thr = 3.0;
    mu2 = 0.4;
    lambda_thr_re = 2.5;
    mu2_re = 3.6;
end

% Ensure sigma_psd resized to Nf is usable by convolving it a bit
if(profile.Nf > 0)
    psd_blur = process_psd_for_nf(sigma_psd, [], profile);
    psd_k = fspecial('gaussian',[1+2*(floor(0.5*size(psd_blur, 1)/profile.Nf(1))), 1], ...
                                 1+2*(floor(0.5*size(psd_blur, 1)/profile.Nf(1))) / 20) * ...
                                 fspecial('gaussian',[1, 1+2*(floor(0.5*size(psd_blur, 2)/profile.Nf(2)))], ...
                                 1+2*(floor(0.5*size(psd_blur, 2)/profile.Nf(2))) / 20);
    
    psd_k = repmat(psd_k, [1, 1, 1+2*(floor(0.5*size(psd_blur, 3)/profile.Nf(3)))]) .* permute(...
         fspecial('gaussian',[1+2*(floor(0.5*size(psd_blur, 3)/profile.Nf(3))), 1], ...
                              1+2*(floor(0.5*size(psd_blur, 3)/profile.Nf(3))) / 20), ...
                                 [2, 3, 1]);
    psd_k = psd_k/sum(psd_k(:));

    psd_blur = convn(padarray(psd_blur, (size(psd_k) - 1) / 2, 'circular'), psd_k, 'valid');
else
    psd_blur = sigma_psd;
    psd_k = 1;
end

% Replace things which had no value previously
if profile.lambda_thr == profile.NO_VALUE; profile.lambda_thr = lambda_thr; end
if profile.mu2 == profile.NO_VALUE; profile.mu2 = mu2; end
if profile.lambda_thr_re == profile.NO_VALUE; profile.lambda_thr_re = lambda_thr_re; end
if profile.mu2_re == profile.NO_VALUE; profile.mu2_re = mu2_re; end

profile.lambda_thr = profile.lambda_thr * profile.filter_strength;
profile.mu2 = profile.mu2 * profile.filter_strength.^2;
profile.lambda_thr_re = profile.lambda_thr_re * profile.filter_strength;
profile.mu2_re = profile.mu2_re * profile.filter_strength.^2;
end



function [t_forward, t_inverse, hadper_trans_single_den, inverse_hadper_trans_single_den, Wwin3D] = ...
    get_transforms(profile, stage_ht)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create transform matrices, etc.
%%%%

if stage_ht
    % get (normalized) forward and inverse transform matrices
    [t_forward_x, t_inverse_x] = get_transf_matrix(profile.N1(1), profile.transform_3D_HT_name, profile.decLevel); 
    [t_forward_y, t_inverse_y] = get_transf_matrix(profile.N1(2), profile.transform_3D_HT_name, profile.decLevel); 
    [t_forward_z, t_inverse_z] = get_transf_matrix(profile.N1(3), profile.transform_3D_HT_name, profile.decLevel); 
else
    % get (normalized) forward and inverse transform matrices
    [t_forward_x, t_inverse_x] = get_transf_matrix(profile.N1_wiener(1), profile.transform_3D_Wiener_name, 0); 
    [t_forward_y, t_inverse_y] = get_transf_matrix(profile.N1_wiener(2), profile.transform_3D_Wiener_name, 0); 
    [t_forward_z, t_inverse_z] = get_transf_matrix(profile.N1_wiener(3), profile.transform_3D_Wiener_name, 0); 

end

t_forward = {permute(single(t_forward_x), [2, 1]), permute(single(t_forward_y), [2, 1]), permute(single(t_forward_z), [2, 1])};
t_inverse = {permute(single(t_inverse_x), [2, 1]), permute(single(t_inverse_y), [2, 1]), permute(single(t_inverse_z), [2, 1])};

if ((strcmp(profile.transform_NL_name, 'haar') == 1) || ...
        (strcmp(profile.transform_NL_name(end-2:end), '1.1') == 1))
    %%% If Haar is used in the NL dimension, then a fast internal transform is used,
    %%% thus no need to generate transform matrices.
    hadper_trans_single_den         = {};
    inverse_hadper_trans_single_den = {};
else
    %%% Create transform matrices. The transforms are later applied by
    %%% matrix-vector multiplication for the 1D case.
    for hpow = 0:ceil(log2(max(profile.N2,profile.N2_wiener)))
        h = 2^hpow;
        [Tfor3rd, Tinv3rd]   = get_transf_matrix(h, profile.transform_NL_name, 0);
        hadper_trans_single_den{h}         = single(Tfor3rd)';
        inverse_hadper_trans_single_den{h} = single(Tinv3rd)';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Kaiser windows used in the aggregation of block-wise estimates
%%%%
% hardcode the window function so that the signal processing toolbox is not needed by default
if profile.beta_wiener==2 && profile.beta==2 && profile.N1_wiener(1) == 8 && profile.N1(1) == 8 && ...
        profile.N1_wiener(2) == 8 && profile.N1(2) == 8 && ...
        profile.N1_wiener(3) == 1 && profile.N1(3) == 1
    Wwin3D = [ 0.1924    0.2989    0.3846    0.4325    0.4325    0.3846    0.2989    0.1924;
        0.2989    0.4642    0.5974    0.6717    0.6717    0.5974    0.4642    0.2989;
        0.3846    0.5974    0.7688    0.8644    0.8644    0.7688    0.5974    0.3846;
        0.4325    0.6717    0.8644    0.9718    0.9718    0.8644    0.6717    0.4325;
        0.4325    0.6717    0.8644    0.9718    0.9718    0.8644    0.6717    0.4325;
        0.3846    0.5974    0.7688    0.8644    0.8644    0.7688    0.5974    0.3846;
        0.2989    0.4642    0.5974    0.6717    0.6717    0.5974    0.4642    0.2989;
        0.1924    0.2989    0.3846    0.4325    0.4325    0.3846    0.2989    0.1924];
elseif profile.beta==2 && all(profile.N1 == [4, 4, 4]) && stage_ht
    kai = [0.438676279837049;0.924313875626079;0.924313875626079;0.438676279837049];
    Wwin3D = repmat(kai * kai', 1, 1, 4) .* permute(kai, [3, 2, 1]);
elseif profile.beta==2 && all(profile.N1_wiener == [5, 5, 5]) && ~stage_ht
    kai = [0.438676279837049;0.834761433401167;1;0.834761433401167;0.438676279837049];
    Wwin3D = repmat(kai * kai', 1, 1, 5) .* permute(kai, [3, 2, 1]);
else 
    if stage_ht
         % Kaiser window used in the aggregation of the HT part
        Wwin3D = kaiser(profile.N1(1), profile.beta) * kaiser(profile.N1(2), profile.beta)';
        Wwin3D = repmat(Wwin3D, 1, 1, profile.N1(3)) .* ...
            permute(kaiser(profile.N1(3), profile.beta), [3, 2, 1]);
    else
        % Kaiser window used in the aggregation of the Wiener filt. part
        Wwin3D = kaiser(profile.N1_wiener(1), profile.beta_wiener) * kaiser(profile.N1_wiener(2), profile.beta_wiener)';
        Wwin3D = repmat(Wwin3D, 1, 1, profile.N1_wiener(3)) .* ...
            permute(kaiser(profile.N1_wiener(3), profile.beta_wiener), [3, 2, 1]);
    end
end

end


function P_rot_f = project_psd_for_param_est(sigma_psd)
% Return the 2-D projection of maximum variance squeezed along a principal
% axis.
% sigma_psd should be a NxNxN matrix.

N=size(sigma_psd,1);
[G2, G1, G3]=meshgrid(1:N,1:N,1:N);

trapz2D=@(G3,G2,G1,sigma_psd) trapz3(G2(:,:,1),trapz3(G3,sigma_psd,3),2);
trapz3D=@(G3,G2,G1,sigma_psd) trapz(G1(:,1,1), trapz2D(G3, G2, G1, sigma_psd)); 

Pn=sigma_psd/trapz3D(G3,G2,G1,sigma_psd);

m2=trapz3D(G3,G2,G1,Pn.*G2);
m1=trapz3D(G3,G2,G1,Pn.*G1);
m3=trapz3D(G3,G2,G1,Pn.*G3);

C=zeros(3);
O1=[2 1 1 1 0 0 1 0 0];
O2=[0 1 0 1 2 1 0 1 0];
O3=[0 0 1 0 0 1 1 1 2];
for jj=[1, 2, 3, 5, 6, 9]
    C(jj)=squeeze(trapz3D(G3,G2,G1,Pn.*(G2-m2).^O2(jj).*(G1-m1).^O1(jj).*(G3-m3).^O3(jj)));
end
C(4) = C(2);
C(7) = C(3);
C(8) = C(6);

% Select highest variance
[U, ~, ~]=svd(C);
P_rot_f = 0;
P_rot_var_f = 0;

for ax_ix=1:3
    v = (U(:, ax_ix));
    P_rot = ifftshift(rotateTo(sigma_psd, v));
    P_rot_proj = squeeze(P_rot(1, :, :));
    P_var = sum(P_rot_proj(:));
    if P_var > P_rot_var_f
        P_rot_f = P_rot_proj;
        P_rot_var_f = P_var;
    end
end
end

function P_rot = rotateTo(sigma_psd, v)
% Rotate sigma_psd to be oriented by v (3x1)
% sigma_psd should be a NxNxN matrix;
% v should be unit length.

r = [1, 0, 0];
if sum(abs(v(:)-r(:))) < 1e-3 || sum(abs(v(:)+r(:))) < 1e-3
    P_rot = sigma_psd;
    return;
end
rotAx = cross(v, r);
rotAx = rotAx / norm(rotAx(:));
rotAngle = acos(min(max(dot(v, r), -1), 1));

rotAxCr = [0, -rotAx(3), rotAx(2); rotAx(3), 0, -rotAx(1); -rotAx(2), rotAx(1), 0];
rotationMatrix = cos(rotAngle) * eye(3) + sin(rotAngle) * rotAxCr + ...
    (1 - cos(rotAngle)) * (rotAx' * rotAx);

N = size(sigma_psd, 1);
N3=3*N;
[G13N, G23N, G33N]=ndgrid((1:N3)-(N3+1)/2, (1:N3)-(N3+1)/2, (1:N3)-(N3+1)/2);
gns = cat(4, G13N, G23N, G33N);
Grot=zeros(N, N, N, 3);
for x=1:N
    for y=1:N
        for z=1:N
        rot = rotationMatrix * squeeze(gns(N+z,N+y,N+x, :));
        Grot(z, y, x, :) = rot;
        end
    end
end
P_rot_handle = griddedInterpolant(G13N,G23N,G33N,repmat(sigma_psd,[3,3,3]),'linear','nearest');
P_rot = P_rot_handle(Grot(:, :, :, 1), Grot(:, :, :, 2), Grot(:, :, :, 3));

end

function I = trapz3(X,Y,dimm)

if dimm==3
    I=sum((Y(:,:,2:end)+Y(:,:,1:end-1))/2.*(X(:,:,2:end)-X(:,:,1:end-1)),3);
elseif dimm==2
    I=sum((Y(:,2:end,:)+Y(:,1:end-1,:))/2.*(X(:,2:end,:)-X(:,1:end-1,:)),2);
else
    I=sum((Y(2:end,:,:)+Y(1:end-1,:,:))/2.*(X(2:end,:,:)-X(1:end-1,:,:)),1);
end

end

function qshifts = get_shift_params(block_size, step_size)
% Select pre-optimized block shifts based on the 1st and 3rd dimension block and step sizes,
% avoiding boundary overlaps.

block_size_1 = min(block_size(1), block_size(2));
block_size_3 = block_size(3);

step_size_1 = min(step_size(1), step_size(2));
step_size_3 = step_size(3);


shift_list = {0, [0, 1], [0, 1, 2, 1], [0, 2, 1], [0, 2, 0, 1], [0, 2], [0, 1, 2], [0, 3, 2, 1], [0, 3], ...
              [0, 2, 3, 1], [0, 1, 1, 0], [0, 3, 1], [0, 2, 1, 3], [0, 4, 3, 1], [0, 3, 4, 1], [0, 4, 2], ...
              [0, 4], [0, 2, 4, 2], [0, 2, 4], [0, 4, 1, 3], [0, 1, 2, 3], [0, 3, 1, 2], [0, 3, 0, 2], ...
              [0, 3, 1, 4], [0, 5, 1, 4], [0, 5], [0, 1, 4, 5], [0, 5, 1], [0, 1, 3], [0, 3, 5, 2], ...
              [0, 5, 2, 3], [0, 1, 3, 4], [0, 5, 4, 1], [0, 2, 4, 1], [0, 3, 1, 5], [0, 6, 2, 4], ...
              [0, 1, 4, 2], [0, 4, 6, 2], [0, 6, 3], [0, 2, 3, 5], [0, 1, 2, 4], [0, 1, 3, 5], ...
              [0, 2, 4, 6], [0, 5, 2], [0, 4, 0, 1], [0, 5, 0, 2], [0, 6, 4, 2], [0, 7, 2, 5], ...
              [0, 6], [0, 1, 4], [0, 5, 3, 2], [0, 2, 5], [0, 5, 7, 2], [0, 7, 2], [0, 7, 1, 6], ...
              [0, 4, 5, 1], [0, 1, 5], [0, 3, 2, 5], [0, 4, 2, 6], [0, 3, 6, 1], [0, 6, 1], ...
              [0, 1, 2, 6], [0, 6, 2], [0, 3, 6], [0, 5, 2, 7], [0, 1, 6, 7], [0, 6, 7, 1], ...
              [0, 2, 7], [0, 7, 6, 1], [0, 1, 7], [0, 3, 6, 3]};

shifts_1 = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 6, 6, 6, 0, 0, 0, 0, 0, 7, 7, 7, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 1, 10, 1, 0, 0, 0, 0, 6, 6, 6, 6, 0, 0, 0, 0, 7, 7, 7, 7, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 6, 6, 6, 0, 0, 0, 0, 0, 5, 5, 5, 0, 0, 0, 0, 0, 5, ...
            5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 1, 10, 10, 0, 0, 0, 0, 6, 6, 6, 6, 0, 0, 0, 0, 5, 5, 5, 5, 0, 0, 0, 0, 5, ...
            5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 1, 10, 10, 10, 0, 0, 0, 6, 6, 6, 6, 6, 0, 0, 0, 5, 5, 5, 5, 5, 0, 0, 0, 5, ...
            5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 10, 1, 0, 0, 0, 0, 0, 3, 3, 3, 0, 0, 0, 0, 0, 20, 20, 20, 0, 0, 0, 0, 0, 22, ...
            5, 5, 0, 0, 0, 0, 0, 15, 15, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 3, 3, 3, 3, 0, 0, 0, 0, 20, 20, 20, 20, 0, 0, 0, 0, 22, ...
            22, 5, 5, 0, 0, 0, 0, 15, 15, 15, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 1, 10, 1, 1, 0, 0, 0, 3, 3, 3, 3, 3, 0, 0, 0, 20, 20, 20, 20, 20, 0, 0, 0, 22, ...
            22, 5, 5, 5, 0, 0, 0, 15, 15, 15, 15, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 3, 3, 3, 3, 3, 3, 0, 0, 20, 20, 20, 20, 20, 20, 0, 0, 22, ...
            22, 22, 5, 5, 5, 0, 0, 15, 15, 15, 15, 15, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 6, 6, 6, 0, 0, 0, 0, 0, 5, 5, 5, 0, 0, 0, 0, 0, 22, ...
            5, 5, 0, 0, 0, 0, 0, 15, 8, 8, 0, 0, 0, 0, 0, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 1, 10, 10, 0, 0, 0, 0, 6, 6, 6, 6, 0, 0, 0, 0, 5, 5, 5, 5, 0, 0, 0, 0, 22, ...
            22, 5, 5, 0, 0, 0, 0, 8, 8, 8, 8, 0, 0, 0, 0, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 1, 10, 10, 10, 0, 0, 0, 6, 6, 6, 6, 6, 0, 0, 0, 5, 5, 5, 5, 5, 0, 0, 0, 22, ...
            22, 5, 5, 5, 0, 0, 0, 8, 8, 8, 8, 8, 0, 0, 0, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 10, 1, 10, 10, 10, 10, 0, 0, 6, 6, 6, 6, 6, 6, 0, 0, 5, 5, 5, 5, 5, 5, 0, 0, 22, ...
            22, 22, 5, 5, 5, 0, 0, 15, 15, 8, 8, 8, 8, 0, 0, 8, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 1, 10, 10, 10, 10, 10, 0, 6, 6, 6, 6, 6, 6, 6, 0, 5, 5, 5, 5, 5, 5, 5, 0, 22, ...
            22, 28, 5, 5, 5, 5, 0, 15, 15, 8, 8, 8, 8, 8, 0, 8, 8, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 0, 0, 0, 0, 0, 6, 6, 6, 0, 0, 0, 0, 0, 9, 9, 9, 0, 0, 0, 0, 0, 44, ...
            44, 44, 0, 0, 0, 0, 0, 8, 8, 8, 0, 0, 0, 0, 0, 45, 5, 5, 0, 0, 0, 0, 0, 46, 46, 46, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 10, 10, 10, 0, 0, 0, 0, 6, 6, 6, 6, 0, 0, 0, 0, 9, 9, 9, 9, 0, 0, 0, 0, 44, ...
            44, 44, 44, 0, 0, 0, 0, 8, 8, 8, 8, 0, 0, 0, 0, 45, 5, 5, 5, 0, 0, 0, 0, 46, 46, 46, 46, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 1, 10, 10, 10, 0, 0, 0, 6, 6, 6, 6, 6, 0, 0, 0, 9, 9, 9, 9, 9, 0, 0, 0, 44, ...
            44, 1, 44, 44, 0, 0, 0, 8, 8, 8, 8, 8, 0, 0, 0, 5, 5, 5, 5, 5, 0, 0, 0, 46, 46, 46, 46, 46, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 10, 10, 10, 0, 0, 3, 3, 6, 3, 3, 6, 0, 0, 9, 9, 9, 9, 9, 9, 0, 0, 44, ...
            44, 44, 44, 44, 44, 0, 0, 8, 8, 8, 8, 8, 8, 0, 0, 45, 45, 5, 5, 5, 5, 0, 0, 46, 46, 46, 46, 46, 46, 0, ...
            0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 10, 10, 10, 10, 10, 10, 10, 0, 3, 3, 3, 3, 3, 3, 6, 0, 9, 9, 9, 9, 9, 9, 9, 0, 44, ...
            44, 44, 6, 44, 44, 44, 0, 8, 8, 8, 8, 8, 8, 8, 0, 45, 45, 5, 5, 5, 5, 5, 0, 46, 46, 46, 46, 46, 46, 46, ...
            0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 1, 10, 10, 10, 10, 10, 10, 10, 3, 6, 3, 6, 3, 3, 3, 6, 9, 9, 9, 9, 9, 9, 9, 9, 44, ...
            44, 44, 44, 44, 44, 44, 44, 8, 8, 8, 8, 8, 8, 8, 8, 45, 45, 5, 5, 5, 5, 5, 5, 46, 46, 46, 46, 46, 46, ...
            46, 46, 71];

shifts_2 = [0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0, 4, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0, 4, 5, 5, 0, 0, 0, 0, 0, 7, 8, 8, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 9, 7, 8, 7, 0, 0, 0, 0, 6, 6, 5, 5, 0, 0, 0, 0, 9, 11, 8, 8, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0, 3, 5, 5, 0, 0, 0, 0, 0, 12, 5, 5, 0, 0, 0, 0, 0, 13, ...
            8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 9, 7, 8, 8, 0, 0, 0, 0, 6, 6, 5, 5, 0, 0, 0, 0, 7, 5, 5, 5, 0, 0, 0, 0, 14, ...
            15, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 7, 7, 8, 8, 8, 0, 0, 0, 15, 15, 16, 16, 16, 0, 0, 0, 17, 18, 5, 5, 5, 0, 0, 0, 19, ...
            15, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0, 4, 5, 5, 0, 0, 0, 0, 0, 21, 8, 8, 0, 0, 0, 0, 0, 23, ...
            8, 8, 0, 0, 0, 0, 0, 24, 16, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 9, 20, 20, 20, 0, 0, 0, 0, 3, 3, 5, 5, 0, 0, 0, 0, 12, 20, 8, 8, 0, 0, 0, 0, 13, ...
            18, 8, 8, 0, 0, 0, 0, 26, 27, 16, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 20, 9, 1, 20, 20, 0, 0, 0, 18, 18, 5, 5, 5, 0, 0, 0, 20, 28, 8, 8, 8, 0, 0, 0, 18, ...
            18, 8, 8, 8, 0, 0, 0, 24, 15, 16, 16, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 9, 20, 20, 20, 20, 20, 0, 0, 15, 18, 18, 5, 5, 5, 0, 0, 29, 30, 30, 25, 25, 25, 0, ...
            0, 31, ...
            23, 18, 8, 8, 8, 0, 0, 32, 24, 27, 16, 16, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0, 12, 5, 5, 0, 0, 0, 0, 0, 33, ...
            5, 5, 0, 0, 0, 0, 0, 34, 16, 16, 0, 0, 0, 0, 0, 35, 25, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 9, 7, 8, 8, 0, 0, 0, 0, 6, 6, 1, 1, 0, 0, 0, 0, 7, 5, 5, 5, 0, 0, 0, 0, 36, ...
            18, 5, 5, 0, 0, 0, 0, 18, 18, 16, 16, 0, 0, 0, 0, 37, 38, 25, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 7, 7, 8, 8, 8, 0, 0, 0, 6, 3, 1, 1, 1, 0, 0, 0, 12, 12, 5, 5, 5, 0, 0, 0, 18, ...
            18, 5, 5, 5, 0, 0, 0, 15, 15, 16, 16, 16, 0, 0, 0, 35, 38, 25, 25, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, ...
            0, 0, 0, 0, 0, 0, 0, 32, 7, 32, 8, 8, 8, 0, 0, 3, 6, 6, 1, 1, 1, 0, 0, 39, 29, 17, 5, 5, 5, 0, 0, 40, ...
            33, 18, 5, 5, 5, 0, 0, 41, 34, 18, 16, 16, 16, 0, 0, 37, 35, 38, 25, 25, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 7, 9, 32, 8, 8, 8, 8, 0, 6, 9, 3, 1, 1, 1, 1, 0, 29, 39, 42, 42, 42, 42, 42, 0, 42, ...
            42, 38, 5, 5, 5, 5, 0, 37, 37, 15, 16, 16, 16, 16, 0, 35, 35, 43, 25, 25, 25, 25, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0, 3, 6, 6, 0, 0, 0, 0, 0, 7, 5, 5, 0, 0, 0, 0, 0, 19, ...
            16, 16, 0, 0, 0, 0, 0, 29, 8, 8, 0, 0, 0, 0, 0, 37, 25, 25, 0, 0, 0, 0, 0, 47, 48, 48, 0, 0, 0, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 9, 7, 8, 8, 0, 0, 0, 0, 6, 6, 6, 6, 0, 0, 0, 0, 20, 11, 5, 5, 0, 0, 0, 0, 49, ...
            49, 16, 16, 0, 0, 0, 0, 50, 51, 8, 8, 0, 0, 0, 0, 35, 51, 25, 25, 0, 0, 0, 0, 52, 53, 48, 48, 0, 0, 0, ...
            0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 7, 9, 8, 8, 8, 0, 0, 0, 15, 15, 15, 15, 15, 0, 0, 0, 12, 11, 5, 5, 5, 0, 0, 0, 49, ...
            6, 1, 16, 16, 0, 0, 0, 29, 39, 8, 8, 8, 0, 0, 0, 51, 43, 25, 25, 25, 0, 0, 0, 54, 53, 48, 48, 48, 0, 0, ...
            0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 55, 26, 32, 8, 8, 8, 0, 0, 27, 27, 15, 56, 56, 15, 0, 0, 39, 29, 51, 5, 5, 5, 0, 0, ...
            26, ...
            55, 55, 16, 16, 16, 0, 0, 39, 57, 8, 8, 8, 8, 0, 0, 58, 42, 51, 25, 25, 25, 0, 0, 52, 47, 54, 48, 48, ...
            48, 0, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 26, 26, 32, 8, 8, 8, 8, 0, 56, 27, 56, 56, 56, 56, 15, 0, 59, 59, 60, 48, 48, 48, ...
            48, 0, 61, ...
            61, 62, 48, 48, 48, 48, 0, 29, 29, 63, 8, 8, 8, 8, 0, 37, 37, 51, 25, 25, 25, 25, 0, 47, 64, 53, 48, 48, ...
            48, 48, 0, 0, ...
            0, 0, 0, 0, 0, 0, 0, 20, 55, 32, 32, 8, 8, 8, 8, 56, 15, 27, 15, 56, 56, 56, 15, 65, 66, 67, 54, 48, 48, ...
            48, 48, 66, ...
            68, 69, 66, 48, 48, 48, 48, 29, 50, 38, 70, 8, 8, 8, 8, 37, 46, 51, 51, 25, 25, 25, 25, 64, 52, 67, 54, ...
            48, 48, 48, 48, 71];

if block_size_3 == 1 && step_size_3 == 1 && 2 < block_size_1 && block_size_1 < 9 && step_size_1 < 8
    shift_list = {0, [0, 1], [0, 2, 1], [0, 2], [0, 1, 2], [0, 3, 2, 1], [0, 3], [0, 1, 1, 0], [0, 1, 2, 3], ...
                  [0, 5], [0, 4], [0, 2, 3, 1], [0, 4, 0, 1], [0, 6]};
    shifts_1 = [0, 1, 3, 0, 0, 0, 0, 0, 0, 1, 3, 6, 0, 0, 0, 0, 0, 1, 3, 3, 6, 0, 0, 0, 0, 1, 3, 6, 6, 9, 0, 0, 0, ...
                1, 1, 3, 3, 10, 9, 0, 0, 1, 4, 3, 10, 6, 9, 13, 14];
    block_size_1 = block_size_1 - 1;
    step_size_1 = step_size_1 - 1;
    ix = (block_size_1 - 2) * 8 + step_size_1;
    sh_1 = shift_list{shifts_1(ix + 1) + 1};
    qshifts = {int32(0), int32(sh_1), int32(0)};
    return;
end

if block_size_1 > 8 || block_size_3 > 8 || step_size_1 > block_size_1 || ...
    step_size_3 > block_size_3 || block_size_3 > block_size_1 || ...
    min(block_size_1, block_size_3) < 3
    qshifts = {int32(0), int32(0), int32(0)};
    return;
end

block_size_1 = block_size_1 - 1;
step_size_1 = step_size_1 - 1;
block_size_3 = block_size_3 - 1;
step_size_3 = step_size_3 - 1;

ix = (block_size_1 - 2) * 6 * 8 * 8 + (block_size_3 - 2) * 8 * 8 + step_size_1 * 8 + step_size_3;
sh_1 = shift_list{shifts_1(ix + 1) + 1};
sh_2 = shift_list{shifts_2(ix + 1) + 1};

qshifts = {int32(sh_2), int32(sh_1), int32(0)};
end