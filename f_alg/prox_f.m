function [X_est, S, C] = prox_f(X, O, Zw, S, C, rhoList, varargin)

    MaxIter = 10;
    Th = 1e-5;

    [I,J,K] = size(X);
    R = size(S,3);

    Ovec = find(O);
    N_ovec = length(Ovec);
    Ovec_inv = find(1-O);

    X_mat = reshape(X,[I*J,K])';
    S_mat = reshape(S,[I*J,R]);
    Z_mat = reshape(Zw,[I*J,R]);
    
    if ~isempty(varargin)
        xi = varargin{1};
    else
        xi_preset = 1e-7;
        xi = xi_preset;
    end
    %% Update unobserved indices
    S_mat(Ovec_inv,:) = Z_mat(Ovec_inv,:);
    S_Omat_temp = S_mat(Ovec,:);

    %% Update observed indices
    X_last = 0;
    for i = 1:MaxIter
        % Update S_r, observed
        temp = ( X_mat(:,Ovec) - C * S_mat(Ovec,:)' );
        for r = 1:R
            rho = rhoList(r);

            temp = temp + C(:,r) * S_mat(Ovec,r)';

            S_r_min = 0;% min(min(S_mat(Ovec,r)));
            S_mat(Ovec,r) = (   temp'*C(:,r) + ...
                rho/2 * Z_mat(Ovec,r)   ) / ...
                ( C(:,r)'*C(:,r) + rho/2 );
            S_mat(Ovec,r) = max(S_mat(Ovec,r),S_r_min);

            C(:,r) = (   temp * S_mat(Ovec,r) / ...
                (   S_mat(Ovec,r)'*S_mat(Ovec,r) + xi   )   );
            C(:,r) = max(C(:,r),0);
            temp = temp - C(:,r) * S_mat(Ovec,r)';
        end

        S = reshape(S_mat,[I,J,R]);
        X_est = reshape(S_mat*C',[I,J,K]);
        dist(i) = norm(X_est(:) - X_last(:))/norm(X_last(:));

        if dist(i) <= Th && i >= 5
            break
        end
        X_last = X_est;
    end

end