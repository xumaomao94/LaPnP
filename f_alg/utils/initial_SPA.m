function [S_est, C_est, R_est] = initial_SPA(X_observed, Omega, R)

[I,J,K] = size(X_observed);
Ovec = find(Omega(:));
Oinv_vec= find(1-Omega(:));

if K == 1 % for matrix decomposition
    S_est = reshape(X_observed,[I*J,1]);
    C_est = 1;
else
    X_mat = reshape(X_observed,I*J,K);
    X_obs_mat = X_mat(Ovec,:); % NUM(O) * K
    
    Normalizer = sum(X_obs_mat,1)';
    X_obs_mat_norm = X_obs_mat * diag( 1./Normalizer );
    
    indices_S = SPA_pnp(X_obs_mat_norm, R);

    R = min(length(indices_S),R);
    R_est = R;
    
    S_obs_mat = X_obs_mat_norm(:,indices_S);
    S_obs_mat = S_obs_mat * diag( Normalizer(indices_S) ); % I * R, or num(O)*R
    
    
        % pseudo_inverse_S = helper.pseudo_inverse(Sm_omega);
        % C = pseudo_inverse_S*Xm_omega;
        
    C_est = (S_obs_mat'*S_obs_mat) \ (S_obs_mat' * X_obs_mat);
    C_est = C_est'; % J * R, or K * R
    S_est = zeros(I*J,R);
    S_est(Ovec,:) = S_obs_mat; % IJ * R
    
    [C_est,Normalizer] = ColumnNormalization(C_est);
        % C_est = C_est * diag(1./Normalizer);
    for rr = 1:R
        if Normalizer(rr) == 0
            S_est(Ovec,rr) = S_est(Ovec,rr);
            S_est(Oinv_vec,rr) = S_est(Oinv_vec,rr);
        else
            S_est(Ovec,rr) = S_est(Ovec,rr) * Normalizer(rr);
            % cnorm = (C_est(:,rr)'*C_temp(:,rr))/(C_est(:,rr)'*C_est(:,rr));
            % S_est(Oinv_vec,rr) = S_est(Oinv_vec,rr) / cnorm;
        end
    end
end
