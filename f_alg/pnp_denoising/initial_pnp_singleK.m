function [S_est, C_est, R_est] = initial_pnp_singleK(X_observed, Omega_aug, R, varargin)
    if isempty(varargin)
        init_method = "interpolation"; %"Inter-Kernel"; % fill the blank using interpolation "TPS";%
    else
        init_method = varargin{1};
    end

    [I,J,K] = size(X_observed);
    if K ~= 1
        error("K must equal 1");
    end
    S_est = zeros(I,J,R);
    
    if size(Omega_aug,3) == 1
        Omega = Omega_aug;
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
        X_interp = interp_random(X_observed.*Omega,Omega,"linear");%"TPS");%
        [U,L,V] = svd(X_interp);
        R_est = R;

        C_est = zeros(1,R);
        for r = 1:R
            S_est(:,:,r) = (U(:,r)*V(:,r)') .* Omega;
            C_est(r) = L(r,r);
        end
        


        [~, S_est, C_est] = prox_f(X_observed, Omega,  zeros([I,J,R]), ...
            S_est, C_est, 0);

        % here remains matrix completion problem
        switch init_method
            case "interpolation" % linear kernel
                for r = 1:R
                    if size(Omega_aug,3) == 2
                        Sr = interp_random(S_est(:,:,r).*realOmega,realOmega,"linear");%"TPS");%
                        Sr(ind_B) = 0;
                        S_est(:,:,r) = Sr;
                    else
                        S_est(:,:,r) = interp_random(S_est(:,:,r),Omega,"linear");%"TPS");%
                    end
                end
            case "TPS" % linear kernel
                for r = 1:R
                    if size(Omega_aug,3) == 2
                        Sr = interp_random(S_est(:,:,r).*realOmega,realOmega,"TPS");%"TPS");%
                        Sr(ind_B) = 0;
                        S_est(:,:,r) = Sr;
                    else
                        S_est(:,:,r) = interp_random(S_est(:,:,r),Omega,"TPS");%"TPS");%
                    end
                    
                end
            case "GaussianKernel"
                for r = 1:R
                    if size(Omega_aug,3) == 2
                        Sr = fit_square(S_est(:,:,r).*realOmega,realOmega);
                        Sr(ind_B) = 0;
                        S_est(:,:,r) = Sr;
                    else
                        S_est(:,:,r) = fit_square(S_est(:,:,r),Omega);
                    end
                end
            case "Inter-Kernel"
                for r = 1:R
                    if size(Omega_aug,3) == 2
                        Sr_K = fit_square(S_est(:,:,r).*realOmega,realOmega);
                        if ~all(Sr_K == 0,'all')
                            Sr_K(ind_B) = 0;
                        end
                        Sr_I = interp_random(S_est(:,:,r),realOmega,"linear");%"TPS");%
                        Sr_I(ind_B) = 0;
                    else
                        Sr_K = fit_square(S_est(:,:,r),Omega);
                        Sr_I = interp_random(S_est(:,:,r),Omega,"linear");%"TPS");%
                    end
                    if K == 1 || all(Sr_K == 0,'all')
                        S_est(:,:,r) = Sr_I;
                    else
                        weight = Sr_K.^2;
                        weight = 0.5 .* weight./max(weight(:));
                        S_est(:,:,r) = weight.*Sr_K + (1-weight).*Sr_I;
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
            case "InterMC"
                X_dist = 2 / (I+1);
                X_grid = linspace( -1+X_dist, 1-X_dist, I)';
                Y_dist = 2 / (J+1);
                Y_grid = linspace( -1+Y_dist, 1-Y_dist, J )';

                % loc=[kron(ones(J,1),X_grid), kron(Y_grid,ones(I,1)) ];
                Ovec = find(Omega);
                for r = 1:R
                    S_r_norm = log10(S_est(:,:,r) + 1e-7);
                    to_add = - min(min(S_r_norm),0);
                    S_r_norm = S_r_norm + to_add;

                    S_est(:,:,r) = slumf_1st_mc_nn_resp(I,X_grid,Y_grid,S_r_norm,0);
                    S_est(:,:,r) = max(10.^(S_est(:,:,r) - to_add) - 1e-7, 0);
                    % S_est(:,:,r) = reshape(10.^(d_est_kt - to_add),I,J);
                end
        end

    end

    % test
    % Ovec = find(Omega);
    % X_est_mat = S_mat*C_est';
    % X_mat = reshape(X_observed,[I*J,K]);
    % sum(abs(X_est_mat(Ovec,:) - X_mat(Ovec,:)),"all") / sum(abs(X_mat(Ovec,:)),'all')

end