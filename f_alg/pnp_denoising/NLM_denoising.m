function [Z_est,denoiser] = NLM_denoising(S,W,rho,lambda,varargin)
    
    R = size(S,3);
    inputNoiseSigma = sqrt(mean(lambda)/rho);  % input noise level

    Z_est = zeros(size(S));
    if ~isempty(varargin)
        denoiser = varargin{1};
    else
        denoiser = cell(1,R);
    end
    for rr = 1:R
        input = S(:,:,rr)+W(:,:,rr);

        input_neg = input(input<0);
        if length(input_neg)>=10 %~isempty(input_neg)
            neg_mid = median(input_neg);
    
            idx_valid = input_neg > neg_mid*5;
            neg_validMin = min( input_neg(idx_valid) );
        else
            neg_validMin = 0;
        end

        input_addition = max(1e-9, 1e-9 - neg_validMin);
        input = max(input + input_addition,1e-9);
        % ----------------------------
        

        Normalizer = max(input(:));
        if Normalizer ~= 0
            input = input/Normalizer;
        end
        [w,h,~]=size(input);
        
        N_window = 5;
        N_p = 2;
        if isempty(varargin)
            [output,denoiser{rr}] = DSG_NLM(input, inputNoiseSigma, N_window, N_p);
        else
            [output,denoiser{rr}] = DSG_NLM(input, inputNoiseSigma, N_window, N_p, denoiser{rr});
        end
        
        Z_est(:,:,rr) = output*Normalizer;
        Z_est(:,:,rr) = max(Z_est(:,:,rr) - input_addition,0);
    end

end