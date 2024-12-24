function [X,t,S,C] = sc_test(alg,X_observed,Omega_aug,R,varargin)


[I,J,K] = size(X_observed);
if ~ismatrix(Omega_aug)
    Omega = Omega_aug(:,:,1);
    BuildingMask = Omega_aug(:,:,2);
    ind_O = find(Omega);
    ind_B = find(BuildingMask);
    ind_realO = setdiff(ind_O,ind_B);
    realOmega = zeros(I,J);
    realOmega(ind_realO) = 1;
else
    Omega = Omega_aug;
    BuildingMask = 0;
    realOmega = Omega;
    ind_realO = find(realOmega);
end


switch alg
    case "interpolation"
        tic
        X = interpolate_slice(X_observed,Omega_aug);
        t = toc;
        S = 0;
        C = 0;
    case "pnp-bm3d"
        if isempty(varargin)
            xi = 1e-7;
            rho = 0.002;
            lambda = rho * 20;
        else
            rho = varargin{1};
            lambda = varargin{2};
            xi = varargin{3};
        end
        [X, t, S, C] = sc_pnp(X_observed, Omega, R, "BM3D", ...
                                'BuildingMask',BuildingMask,...
                                'rho',rho,...
                                'lambda',lambda,....
                                'xi',xi);
    case "pnp-dpir"
        if isempty(varargin)
            xi = 1e-7;
            rho = 0.01;
            lambda = 2.5e-4 * rho;
        else
            rho = varargin{1};
            lambda = varargin{2};
            xi = varargin{3};
        end
        [X, t, S, C] = sc_pnp(X_observed, Omega, R, "DPIRNet", ...
                                'BuildingMask',BuildingMask,...
                                'rho',rho,...
                                'lambda',lambda,...
                                'xi',xi);
    case "pnp-nlm"
        if isempty(varargin)
            xi = 1e-7;
            rho = 0.002;
            lambda = rho * 0.01;
        else
            rho = varargin{1};
            lambda = varargin{2};
            xi = varargin{3};
        end
        [X, t, S, C] = sc_pnp(X_observed, Omega, R, "NLM", ...
                                'BuildingMask',BuildingMask,...
                                'rho',rho,...
                                'xi',xi,...
                                'lambda',lambda);
    case "nasdac"
        Ovec = find(Omega);
        Ov = false(1,numel(Omega));
        Ov(Ovec) = true;
        [X, t, S, C] = nasdac(X_observed, Ov, R, 0);
    case "nasdacRobust" % fill the boundary of X_observed and Omega to make it with size [51,51,k]
        target_size = [51, 51, size(X_observed, 3)];
        original_size = size(X_observed);
        
        % Initialize the figure with the original data
        X_extended = X_observed;
        Omega_extended = Omega;
        
        % Duplicate the original image and Omega to all four sides with mirror images
        while size(X_extended, 1) < target_size(1) || size(X_extended, 2) < target_size(2)
            % Mirror the image horizontally
            X_extended = [fliplr(X_extended), X_extended, fliplr(X_extended)];
            Omega_extended = [fliplr(Omega_extended), Omega_extended, fliplr(Omega_extended)];
            % Mirror the image vertically
            X_extended = [flipud(X_extended); X_extended; flipud(X_extended)];
            Omega_extended = [flipud(Omega_extended); Omega_extended; flipud(Omega_extended)];
        end
        
        % Crop the middle [51, 51, k] section
        start_idx = floor((size(X_extended) - target_size) / 2) + 1;
        end_idx = start_idx + target_size - 1;
        X_observed = X_extended(start_idx(1):end_idx(1), start_idx(2):end_idx(2), :);
        Omega = Omega_extended(start_idx(1):end_idx(1), start_idx(2):end_idx(2), :);

        % Normalization to maximum 0.1
        max_val = max(X_observed(:));
        if max_val ~= 0
            X_observed = X_observed / max_val * 0.1;
        end

        % nasdac
        Ovec = find(Omega);
        Ov = false(1,numel(Omega));
        Ov(Ovec) = true;
        [X, t, S, C] = nasdac(X_observed, Ov, R, 0);

        % Recover the original size from X
        start_idx = floor((size(X) - original_size) / 2) + 1;
        end_idx = start_idx + original_size - 1;
        X = X(start_idx(1):end_idx(1), start_idx(2):end_idx(2), :);
        X = X * max_val / 0.1;
        X = X.*(1-BuildingMask);
        S = S(start_idx(1):end_idx(1), start_idx(2):end_idx(2), :);
    case "dowjons"
        Ovec = find(Omega);
        Ov = false(1,numel(Omega));
        Ov(Ovec) = true;
        [X, t, S, C] = dowjons(X_observed, Ov, R, 0);
    case "dowjonsRobust"
        target_size = [51, 51, size(X_observed, 3)];
        original_size = size(X_observed);
        
        % Initialize the figure with the original data
        X_extended = X_observed;
        Omega_extended = Omega;
        
        % Duplicate the original image and Omega to all four sides with mirror images
        while size(X_extended, 1) < target_size(1) || size(X_extended, 2) < target_size(2)
            % Mirror the image horizontally
            X_extended = [fliplr(X_extended), X_extended, fliplr(X_extended)];
            Omega_extended = [fliplr(Omega_extended), Omega_extended, fliplr(Omega_extended)];
            % Mirror the image vertically
            X_extended = [flipud(X_extended); X_extended; flipud(X_extended)];
            Omega_extended = [flipud(Omega_extended); Omega_extended; flipud(Omega_extended)];
        end
        
        % Crop the middle [51, 51, k] section
        start_idx = floor((size(X_extended) - target_size) / 2) + 1;
        end_idx = start_idx + target_size - 1;
        X_observed = X_extended(start_idx(1):end_idx(1), start_idx(2):end_idx(2), :);
        Omega = Omega_extended(start_idx(1):end_idx(1), start_idx(2):end_idx(2), :);

        % Normalization to maximum 0.1
        max_val = max(X_observed(:));
        if max_val ~= 0
            X_observed = X_observed / max_val * 0.1;
        end

        % dowjons
        Ovec = find(Omega);
        Ov = false(1,numel(Omega));
        Ov(Ovec) = true;
        [X, t, S, C] = dowjons(X_observed, Ov, R, 0);

        % Recover the original size from X
        start_idx = floor((size(X) - original_size) / 2) + 1;
        end_idx = start_idx + original_size - 1;
        X = X(start_idx(1):end_idx(1), start_idx(2):end_idx(2), :);
        X = X * max_val / 0.1;
        X = X.*(1-BuildingMask);
        S = S(start_idx(1):end_idx(1), start_idx(2):end_idx(2), :);
    case "blockterm"
        lambda = 0.001*ones(1,3);%[0.2 0.2 0.2];
        L = 4;
        [X, t, A, B, C, S] = ...
            blockterm(X_observed, Omega, R, L, lambda, X_observed, "fiber");%
    case "TPS"
        tic;
        X = zeros(I,J,K);
        for k = 1:K
            X(:,:,k) = interp_random(X_observed(:,:,k),realOmega,"TPS");%"TPS");%
        end
        S = 0;
        C = 0;
        t = toc;
end



end