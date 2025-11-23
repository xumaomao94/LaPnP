function [y_est, blocks] = BM3D_mac_adapt(z, sigma_psd, profile, stage_arg, blockmatches)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  BM3D is an algorithm for attenuation of additive spatially correlated
%  stationary (aka colored) Gaussian noise in grayscale and multichannel images.
%
%
%  FUNCTION INTERFACE:
%
%  y_est = BM3D(z, sigma_psd, profile)
%
%  INPUT ARGUMENTS:
%
%  -- required --
%
%         'z' : noisy image (M x N or M x N x C double array, intensities in range [0,1])
%               For multichannel images, block matching is performed on the first channel.
%  'sigma_psd' : noise power spectral density (M x N double nonnegative array)
%               OR
%               noise STD
%               OR
%               either of these in multichannel:
%               (M x N x C PSDs or 1 x C STDs)
%
% -- optional --
%
%   'profile' : 'np' --> Normal Profile (default)
%               'refilter' --> Apply refiltering
%               OR
%               a BM3DProfile object specifying the parameters
%               some other premade profiles also included from the previous versions
%               in BM3DProfile.m
%
%   'stage_arg' : Determines whether to perform hard-thresholding or wiener filtering.
%                 either BM3DProfile.HARD_THRESHOLDING, BM3DProfile.ALL_STAGES or an estimate
%                  of the noise-free image.
%                    - BM3DProfile.ALL_STAGES: Perform both.
%                    - BM3DProfile.HARD_THRESHOLDING: Perform hard-thresholding only.
%                    - ndarray, size of z: Perform Wiener Filtering with stage_arg as pilot.
%
%   'blockmatches' : Tuple {HT, Wiener}, with either value either:
%                      - false : Do not save blockmatches for phase
%                      (default)
%                      - true : Save blockmatches for phase
%                      - Pre-computed block-matching array returned by a
%                      previous call with [true]
%  OUTPUT:
%      'y_est'  denoised image  (M x N double array)
%      'y_est', {'blocks_ht', 'blocks_wie'} denoised image, plus HT and
%          Wiener blockmatches, if any storeBM values are set to True
%          (or [0] for missing block array, if only one calculated)
%
%
%  BASIC SIMULATION EXAMPLES:
%
%     Case 1)
%
%      % Read a grayscale noise-free image
%
%      y=im2double(imread('cameraman.tif'));
%
%      % Generate noisy observations corrupted by additive colored random noise
%        generated as convution of AWGN against with kernel 'k'
%
%      k=[-1;2;-1]*[1 4 1]/100;   % e.g., a diagonal kernel
%      z=y+imfilter(randn(size(y)),k(end:-1:1,end:-1:1),'circular');
%
%      % define 'sigma_psd' from the kernel 'k'
%
%      sigma_psd=abs(fft2(k,size(z,1),size(z,2))).^2*numel(z);
%
%      % Denoise 'z'
%      y_est = BM3D(z, sigma_psd);
%
%
%     Case 2)
%
%      % Read a grayscale noise-free image
%
%      y=im2double(imread('cameraman.tif'));
%
%      % Generate noisy observations corrupted by additive colored random noise
%      % generated as convution of AWGN against with kernel 'k'
%      [x2, x1]=meshgrid(ceil(-size(y,2)/2):ceil(size(y,2)/2)-1,ceil(-size(y,1)/2):ceil(size(y,1)/2)-1)
%      sigma_psd=ifftshift(exp(-((x1/size(y,1)).^2+(x2/size(y,2)).^2)*10))*numel(y)/100;
%      z=y+real(ifft2(fft2(randn(size(y))).*sqrt(sigma_psd)/sqrt(numel(y))));
%
%      % Denoise 'z'
%      y_est = BM3D(z, sigma_psd);
%
%     Case 3) If 'sigma_psd' is a singleton, this value is taken as sigma and
%             it is assumed that the noise is white variance sigma^2.
%
%      % Read a grayscale noise-free image
%
%      y=im2double(imread('cameraman.tif'));
%
%      % Generate noisy observations corrupted by additive white Gaussian noise with variance sigma^2
%      sigma=0.1;
%      z=y+sigma*randn(size(y));
%
%      y_est = BM3D(z, sigma);
%
%      % or, equivalently,
%      sigma_psd = ones(size(z))*sigma^2*numel(z)
%      y_est = BM3D(z, sigma_psd)
%
%
%      Case 4)   MULTICHANNEL PROCESSING
%
%      y_est = BM3D(cat(3, z1, z2, z3), sigma_psd, 'np'); 
%
%      Multiple PSDs are optionally handled in the same way.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2006-2019 Tampere University.
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
    profile = BM3DProfile(profile);
end

if ~exist('stage_arg','var')
    stage_arg = profile.ALL_STAGES;  % By default, do both HT and Wie
end

if min(size(z, 1), size(z, 2)) < profile.N1 || min(size(z, 1), size(z, 2)) < profile.N1_wiener
    disp('Error: Image cannot be smaller than block size!')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pro = convertToBM4DProfile(profile);

channel_count = size(z, 3);
dont_do_bm = ~exist('blockmatches','var') || numel(blockmatches) ~= 2 || ...
                (blockmatches{1} == false && blockmatches{2} == false);

if channel_count > 1  % Channel dimension should be 4
    z = permute(z, [1, 2, 4, 3]);
    if ndims(sigma_psd) == 3
        sigma_psd = permute(sigma_psd, [1, 2, 4, 3]);
    end

    if ~dont_do_bm
        disp('Warning: block-match data supplied with multichannel BM3D will be discarded. Call the function separately for each channel.');
    end
    [y_est] = BM4D_multichannel(z, sigma_psd, pro, stage_arg);
    y_est = squeeze(y_est);
    return
end

if dont_do_bm
    [y_est] = BM4D(z, sigma_psd, pro, stage_arg);
else
    [y_est, blocks] = BM4D(z, sigma_psd, pro, stage_arg, blockmatches);
end

end

function pro = convertToBM4DProfile(pro_in)
    pro = BM4DProfile('BM3D');

    pro.filter_strength = pro_in.filter_strength;

    pro.print_info = pro_in.print_info;

    pro.transform_3D_HT_name = pro_in.transform_2D_HT_name;
    pro.transform_3D_Wiener_name = pro_in.transform_2D_Wiener_name;
    pro.transform_NL_name = pro_in.transform_3rd_dim_name;

    pro.Nf = [pro_in.Nf, pro_in.Nf, 1];
    pro.Kin = pro_in.Kin;

    pro.denoise_residual = pro_in.denoise_residual;
    pro.residual_thr = pro_in.residual_thr;
    pro.max_pad_size = pro_in.max_pad_size;

    pro.gamma = pro_in.gamma;

    pro.N1 = [pro_in.N1, pro_in.N1, 1];
    pro.Nstep = [pro_in.Nstep, pro_in.Nstep, 1];

    pro.N2 = pro_in.N2;
    pro.Ns = [floor(pro_in.Ns / 2), floor(pro_in.Ns / 2), 0];
    pro.tau_match = pro_in.tau_match * pro_in.N1^2 / 255^2;

    pro.lambda_thr = pro_in.lambda_thr3D;
    pro.mu2 = pro_in.mu2;

    pro.lambda_thr_re = pro_in.lambda_thr3D_re;
    pro.mu2_re = pro_in.mu2_re;
    pro.beta = pro_in.beta;

    pro.N1_wiener = [pro_in.N1_wiener, pro_in.N1_wiener, 1];
    pro.Nstep_wiener = [pro_in.Nstep_wiener, pro_in.Nstep_wiener, 1];

    pro.N2_wiener = pro_in.N2_wiener;
    pro.Ns_wiener = [floor(pro_in.Ns_wiener / 2), floor(pro_in.Ns_wiener / 2), 0];
    pro.tau_match_wiener = pro_in.tau_match_wiener * pro_in.N1_wiener^2 / 255^2;
    pro.beta_wiener = pro_in.beta_wiener;
    pro.decLevel = pro_in.decLevel;

    pro.set_sharpen(pro_in.sharpen_alpha);
    pro.num_threads = pro_in.num_threads;

end