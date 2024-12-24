function ssim_mean = ssim_score_4sc(X,Y)
    
    if size(X)~=size(Y)
        error('The compared tensor should be with the same size')
    end

    X = reshape(X,[size(X,1),size(X,2),numel(X)/size(X,1)/size(X,2)]);
    [I,J,K] = size(X);
    Y = reshape(Y,[I,J,K]);
    
    X = 10*log10(X+1e-12);
    Y = 10*log10(Y+1e-12);
    % X = 10*log10(X+1e-10);
    % Y = 10*log10(Y+1e-10);

    min_xy = min([Y(:)]);
    X = X + min_xy;
    Y = Y + min_xy;

    max_xy = max([Y(:)]);

    X = X./max_xy.*255;
    Y = Y./max_xy.*255;

    ssim_K = zeros(K,1);

    for k = 1:K
        ssim_K(k) = ssim_index(X(:,:,k), Y(:,:,k));
    end

    ssim_mean = mean(ssim_K);
% figure; plot(ssim_K);
end