function [Z_est,Normalizer] = BM3D_denoising(S,W,rho,lambda,varargin)
    
    R = size(S,3);
    inputNoiseSigma = sqrt(mean(lambda)/rho);  % input noise level

    Z_est = zeros(size(S));
    for rr = 1:R
        input = S(:,:,rr)+W(:,:,rr);
        input_neg = input(input<0);
        if length(input_neg)>=10 
            neg_mid = median(input_neg);
            idx_valid = input_neg > neg_mid*5;
            neg_validMin = max( min( input_neg(idx_valid) ), 4e-9 );
        else
            neg_validMin = 0;
        end
        logadd = max(1e-9, 1e-9 - neg_validMin);
        input = 10*log10(max(input+logadd,1e-9));

        if min(input(:))<0
            input_addition = - min(input(:));
        else
            input_addition = 0;
        end
        input = input + input_addition;
        if isempty(varargin)
            Normalizer = max(input(:));
        else
            Normalizer = varargin{1};
        end
        input = input/Normalizer;
        [w,h,~]=size(input);
        
        [~,output] = BM3D(1, input, inputNoiseSigma);
        
        Z_est(:,:,rr) = output*Normalizer;
        Z_est(:,:,rr) = Z_est(:,:,rr) - input_addition;

        Z_est(:,:,rr) = 10.^(Z_est(:,:,rr)./10);
        Z_est(:,:,rr) = max(Z_est(:,:,rr) - logadd,0);
    end

end