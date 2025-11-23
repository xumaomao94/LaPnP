% BM4DProfile: Specify parameters for BM4D
classdef BM4DProfile
    properties (Constant)
        NO_VALUE = -1;
        HARD_THRESHOLDING = int8(1);
        WIENER_FILTERING = int8(2);
        ALL_STAGES = BM4DProfile.HARD_THRESHOLDING + BM4DProfile.WIENER_FILTERING;
    end
    properties
        
        filter_strength = 1;
        
         %%%% Select transforms ('dct', 'dst', 'hadamard', or anything that is listed by 'help wfilters'):
        transform_3D_HT_name     = 'bior1.5';%'dct';% %% transform used for the HT filt. of size N1 x N1
        transform_3D_Wiener_name = 'dct';     %% transform used for the Wiener filt. of size N1_wiener x N1_wiener
        transform_NL_name   = 'haar';%'dct';%    %% transform used in the 3-rd dim, the same for HT and Wiener filt.

        % -- Exact variances for correlated noise: --

        % Variance calculation parameters
        Nf = [16, 16, 16];  % domain size for FFT computations
        Kin = 4;  % how many layers of var3D to calculate

        denoise_residual = false; % Perform residual thresholding and re-denoising
        residual_thr = 3; % Threshold for the residual HT (times sqrt(PSD))
        max_pad_size = [-1, -1, -1]; % Maximum required pad size (= half of the kernel size), or [-1, -1] -> use image size

        % Block matching
        gamma = 3.0;
        
        print_info = false;
        
        %%%% Hard-thresholding (HT) parameters:
        N1                  = [4, 4, 4];%   %% N1 x N1 is the block size used for the hard-thresholding (HT) filtering
        Nstep               = [3, 3, 3];   %% sliding step to process every next reference block
        N2                  = 16;  %% maximum number of similar blocks (maximum size of the 3rd dimension of a 3D array)
        Ns                  = [7, 7, 7];  %% half length of the side of the search neighborhood for full-search block-matching (BM), must be odd
        tau_match           = 2.9527;%% threshold for the block-distance (d-distance)
        beta                = 2.0; %% parameter of the 2D Kaiser window used in the reconstruction
        
        
        
        % threshold parameter for the hard-thresholding in 3D transform domain
        % NO_VALUE in lambda or lamdba_wiener means automatic selection.
        
        lambda_thr        = BM4DProfile.NO_VALUE; %2.7;
        % Refiltering
        lambda_thr_re       = BM4DProfile.NO_VALUE;
        
        %%%% Wiener filtering parameters:
        N1_wiener           = [5, 5, 5];%4;%
        Nstep_wiener        = [3, 3, 3];
        N2_wiener           = 32;%8;%
        Ns_wiener           = [7, 7, 7];
        tau_match_wiener    = 0.7689;
        beta_wiener         = 2.0;
        
        mu2       = BM4DProfile.NO_VALUE %1.0;
        mu2_re      = BM4DProfile.NO_VALUE;
        
        splitComp = [false, false, false];
        
        decLevel = 0;        %% dec. levels of the dyadic wavelet 2D transform for blocks (0 means full decomposition, higher values decrease the dec. number)
        % Unused
        adjust_complex_params = false;
        % Sharpening parameters, 1=no sharpen, >1, sharpen
        % use set_sharpen to set both. Use HARD_THRESHOLDING stage
        % when sharpening.
        sharpen_alpha = 1;
        sharpen_alpha_3d = 1;
        num_threads = 0; % 0=Automatic. Multithreaded results have slight
        % variations in output due to changing order of aggregation, set
        % num_threads=1 for deterministic output
    end
    methods

        function p = set_sharpen(p, alpha)
            p.sharpen_alpha = alpha;
            p.sharpen_alpha_3d = alpha;
        end
        
        function pro = BM4DProfile(profileName)
            
            if ~exist('profileName','var')
                profileName = 'np'; %% default profile
            end
            
            if strcmp(profileName, 'np') == 1
            elseif strcmp(profileName, 'refilter') == 1
                pro.denoise_residual = true;
                
            elseif strcmp(profileName, 'lc') == 1
                
                pro.Ns                  = [4, 4, 4];
                pro.N2_wiener           = 16;
                pro.Ns_wiener           = [4, 4, 4];
                          
                
            elseif strcmp(profileName, '8x8') == 1
                pro.N1                  = [8, 8, 1];
                pro.Nstep               = [3, 3, 1]; 
                pro.N1_wiener           = [8, 8, 1];
                pro.Nstep_wiener        = [3, 3, 1];   
                
            elseif strcmp(profileName, '8x8refilter') == 1
                pro.N1                  = [8, 8, 1];
                pro.Nstep               = [3, 3, 1]; 
                pro.N1_wiener           = [8, 8, 1];
                pro.Nstep_wiener        = [3, 3, 1];  
                pro.denoise_residual    = true;  
                
            elseif strcmp(profileName, 'complex') == 1
                pro.transform_3D_HT_name     = 'fft';
                pro.transform_3D_Wiener_name = 'fft';
                pro.adjust_complex_params    = true;
                
            elseif strcmp(profileName, 'BM3D') == 1               
                pro.N1                  = [8, 8, 1];
                pro.Nstep               = [3, 3, 1]; 
                pro.N1_wiener           = [8, 8, 1];
                pro.Nstep_wiener        = [3, 3, 1]; 
                pro.Nf                  = [32, 32, 1];
                pro.Ns                  = [19, 19, 1];
                pro.Ns_wiener           = [19, 19, 1];
                
            elseif strcmp(profileName, 'cBM3D') == 1               
                pro.N1                  = [8, 8, 1];
                pro.Nstep               = [3, 3, 1]; 
                pro.N1_wiener           = [8, 8, 1];
                pro.Nstep_wiener        = [3, 3, 1]; 
                pro.Nf                  = [32, 32, 1];
                pro.Ns                  = [19, 19, 1];
                pro.Ns_wiener           = [19, 19, 1];
                pro.transform_3D_HT_name     = 'fft';
                pro.transform_3D_Wiener_name = 'fft';
                pro.adjust_complex_params    = true;
            else
                disp('Error: profile not found! Returning default profile.')
            end
            
        end
    end
end




