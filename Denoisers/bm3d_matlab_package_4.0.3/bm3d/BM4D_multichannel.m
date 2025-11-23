function y_est = BM4D_multichannel(z, sigma_psd, profile, stage_arg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  BM4D is an algorithm for attenuation of additive spatially correlated
%  stationary (aka colored) Gaussian noise in volumetric data.
%  BM4D_multichannel processes 4-D data as multichannel 3-D data,  
%  performing blockmatching only on the first channel.
%  Channel dimension should be the last (such that z = cat(4, ch1, ch2, ...))  
%
%  FUNCTION INTERFACE:
%
%  y_est = BM4D_multichannel(z, sigma_psd, profile)
%
%  INPUT ARGUMENTS:
%
%  -- required --
%
%         'z' : noisy image (4-D array)
%  'sigma_psd' : One noise power spectral density (3-D nonnegative array)
%               OR
%               One noise STD
%               OR one PSD for each channel (4-D array)
%               OR noise STD for each channel (1-D array)
%           
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
        profile = 'np';
    end

    if ~exist('stage_arg','var')
        stage_arg = BM4DProfile.ALL_STAGES;
    end
    
    multiple_sigma = (ndims(sigma_psd) == 1 && numel(sigma_psd) > 1) || ndims(sigma_psd) == 4;

    function sig = get_sigma(idx)
        if ~multiple_sigma
            sig = sigma_psd;
        elseif ndims(sigma_psd) == 1
            sig = sigma_psd(idx);
        else
            sig = sigma_psd(:, :, :, idx);
        end
    end

    y_est = zeros(size(z), 'single');

    [y_est(:, :, :, 1), match_arrs] = BM4D(z(:, :, :, 1), get_sigma(1), profile, stage_arg, {true, true});

    for i = 2:size(z, 4)
        y_est(:, :, :, i) = BM4D(z(:, :, :, i), get_sigma(i), profile, stage_arg, match_arrs);
    end

end

