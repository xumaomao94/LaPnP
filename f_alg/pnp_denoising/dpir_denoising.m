function Z_est = dpir_denoising(S,W,rho,lambda,BuildingMask)
% Set BuildingMask as 0, or zeros(I,J), if there is no

    if all(BuildingMask==0)
        idx_Building = [];
    else
        idx_Building = find(BuildingMask);
    end
    [I_s,J_s,R] = size(S);
    % R = size(S,3);

    inputNoiseSigma = 255*sqrt(lambda/rho);  % input noise level
    
    Z_est = zeros(size(S));
    for rr = 1:R
        input = S(:,:,rr)+W(:,:,rr);

        % ---------------------log input
        input_neg = input(input<0);
        if length(input_neg)>=10 %~isempty(input_neg)
            neg_mid = median(input_neg);
            idx_valid = input_neg > neg_mid*2;
            neg_validMin = max( min( input_neg(idx_valid) ), -4e-9 );
        else
            neg_validMin = 0;
        end
        logadd = max(1e-9, 1e-9 - neg_validMin); % no more than 5e-9
        input = 10*log10(max(input+logadd,1e-9));
        % ----------------------------

        if min(input(:))<0
            input_addition = 0 - min(input(:));
        else
            input_addition = 0;
        end

        input = input + input_addition;
        Normalizer = max(input(:));
        input = input./Normalizer;

        [w,h,~]=size(input);
        
        input = im2single(input);
        if mod(w,2)==1
            input = cat(1,input, input(end,:)) ;
        end
        if mod(h,2)==1
            input = cat(2,input, input(:,end)) ;
        end

        [I,J] = size(input);
        % only for [64,64], [128,128], [256,256], [512,512]
        if I <= 64
            I_max = 64;
        elseif I <= 128
            I_max = 128;
        elseif I <= 256
            I_max = 256;
        elseif I <= 512
            I_max = 512;
        else
            error('size no larger than 512')
        end
        if J <= 64
            J_max = 64;
        elseif I <= 128
            J_max = 128;
        elseif I <= 256
            J_max = 256;
        elseif I <= 512
            J_max = 512;
        else
            error('size no larger than 512')
        end

        I_b = ( I_max - I )/2;
        J_b = ( J_max - J )/2;
        input_ex = zeros(I_max,J_max);
        input_ex(I_b+1:I_b+I, J_b+1:J_b+J) = input;
        input_ex(I_b:-1:1,:) = input_ex(I_b+1:2*I_b,:);
        input_ex(I_b+I+1:end,:) = input_ex(I_b+I:-1:I+1,:);
        input_ex(:,J_b:-1:1) = input_ex(:,J_b+1:2*J_b);
        input_ex(:,J_b+J+1:end) = input_ex(:,J_b+J:-1:J+1);
        
        
        input_py = py.numpy.array(input_ex);
        Sigma_py = py.numpy.array(inputNoiseSigma);

        output = py.dpir_denoising_interface.dpir(input_py,inputNoiseSigma); % matconvnet default
        output = single(output)./255;

    
        output = output(I_b+1:I_b+I,J_b+1:J_b+J);


        if mod(w,2)==1
            output = output(1:end-1,:);
            input  = input(1:end-1,:);
        end
        if mod(h,2)==1
            output = output(:,1:end-1);
            input  = input(:,1:end-1);
        end
        
        output = min(output,1);
        output = max(output,0);
        Z_est(:,:,rr) = single(output)*Normalizer;
        
        Z_est(:,:,rr) = Z_est(:,:,rr) - input_addition;

        % ---------------log input
        Z_est(:,:,rr) = 10.^(Z_est(:,:,rr)./10);
        Z_est(:,:,rr) = max(Z_est(:,:,rr) - logadd,0);
        % ----------------------
    end
end