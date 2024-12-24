function [S_est, C_est, R_est] = initial_pnp(X_observed, Omega_aug, R, rho, xi, varargin)
    if isempty(varargin)
        init_method = "interpolation"; %"Inter-Kernel"; % fill the blank using interpolation "TPS";%
    else
        init_method = varargin{1};
    end

    [I,J,K] = size(X_observed);
    if K == 1
        R = 1;
    end
    S_est = zeros(I,J,R);
    
    if size(Omega_aug,3) == 1
        Omega = Omega_aug;
        realOmega = Omega;
        ind_realO = find(realOmega);
    elseif size(Omega_aug,3) == 2
        Omega = Omega_aug(:,:,1);
        BuildingMask = Omega_aug(:,:,2);
        ind_O = find(Omega);
        ind_B = find(BuildingMask);
        ind_realO = setdiff(ind_O,ind_B);
        realOmega = zeros(I,J);
        realOmega(ind_realO) = 1;
    end

    if ismatrix(Omega) % fiber sampling, then use SPA to recover C
        if K >= R && K ~= 1
            [S_mat, C_est, R_est] = initial_SPA(X_observed, Omega, R);
            R = R_est;
            [C_est,Normalizer] = ColumnNormalization(C_est);
            S_mat = S_mat * diag(Normalizer);
        elseif K == 1
            [S_mat, C_est] = initial_SPA(X_observed, Omega, K);
            R_est = 1;
            R = 1;
            [C_est,Normalizer] = ColumnNormalization(C_est);
            S_mat = S_mat * diag(Normalizer);
        end
        S_est = reshape(S_mat,[I,J,R]);

        [~, S_est, C_est] = prox_f(X_observed, realOmega,  zeros([I,J,R]), ...
            S_est, C_est, rho, xi);
        [C_est,Normalizer] = ColumnNormalization(C_est);
        S_mat = S_mat * diag(Normalizer);

        % completion
        switch init_method
            case "random"
                for r = 1:R
                    Sr = S_est(:,:,r);
                    Smean = mean(Sr(ind_realO));
                    Svar = var(Sr(ind_realO));
                    S_est(:,:,r) = normrnd(Smean,Svar,[I,J]);
                end
            case "interpolation" % linear kernel
                for r = 1:R
                    if size(Omega_aug,3) == 2
                        Sr = interp_random(S_est(:,:,r).*realOmega,realOmega,"linear");
                        S_est(:,:,r) = Sr;
                    else
                        S_est(:,:,r) = interp_random(S_est(:,:,r),Omega,"linear");
                    end
                end
            case "TPS" % linear kernel
                for r = 1:R
                    if size(Omega_aug,3) == 2
                        Sr = interp_random(S_est(:,:,r).*realOmega,realOmega,"TPS");
                        Sr(ind_B) = 0;
                        S_est(:,:,r) = Sr;
                    else
                        S_est(:,:,r) = interp_random(S_est(:,:,r),Omega,"TPS");
                    end
                    
                end
            case "krig"
                options.polytrend=2;
                V='1 Sph(0.25)'; % Select variogram model

                X_dist = 2 / (I+1);
                X_grid = linspace( -1+X_dist, 1-X_dist, I)';
                Y_dist = 2 / (J+1);
                Y_grid = linspace( -1+Y_dist, 1-Y_dist, J )';

                loc=[kron(ones(J,1),X_grid), kron(Y_grid,ones(I,1)) ];
                Ovec = find(Omega);
                for r = 1:R
                    S_r_norm = log10(S_mat(Ovec,r) + 1e-7);
                    to_add = - min(min(S_r_norm),0);
                    S_r_norm = S_r_norm + to_add;
                    % to_devide = max(S_r_norm);
                    % S_r_norm = S_r_norm./to_devide;

                    [d_est_kt,~]=krig(loc(Ovec,:),S_r_norm,loc,V,options);
                    S_est(:,:,r) = reshape(max(10.^(d_est_kt - to_add)-1e-7,0),I,J);
                end
        end

    end

end