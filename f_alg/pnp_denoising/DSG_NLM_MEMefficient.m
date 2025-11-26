function output = DSG_NLM_MEMefficient(input, sqrNoise, N_w, N_p)

    Npower      = sqrNoise^2;
    patch_size  = 2*N_p + 1;
    window_size = 2*N_w + 1;

    [I,J] = size(input);  % only support 2D currently
    N     = I*J;

    % ---------------------------------------------------------------------
    % Pad boundary by replication
    % ---------------------------------------------------------------------
    input_fillupBoundary = zeros(I+2*N_p,J+2*N_p);
    input_fillupBoundary(N_p+1:N_p+I, N_p+1:N_p+J) = input;

    % left/right
    input_fillupBoundary(N_p+1:N_p+I, 1:N_p) = repmat(input(:,1),1,N_p);
    input_fillupBoundary(N_p+1:N_p+I, N_p+J+1:2*N_p+J) = repmat(input(:,end),1,N_p);

    % top/bottom
    input_fillupBoundary(1:N_p,:) = repmat(input_fillupBoundary(N_p+1,:),N_p,1);
    input_fillupBoundary(N_p+I+1:2*N_p+I,:) = repmat(input_fillupBoundary(N_p+I,:),N_p,1);

    % log-domain
    input_fillupBoundary = log10(1e-15 + input_fillupBoundary);

    %% Denoising
    output = zeros(I,J);
    denom = 2 * Npower * (patch_size^2);

    for i = 1:I
        for j = 1:J
            % patch center (in padded coords)
            pi = i;
            pj = j;

            patch_i = input_fillupBoundary(pi:pi+2*N_p, pj:pj+2*N_p);

            % search window limits in original coordinates
            i_min = max(i - N_w, 1);
            i_max = min(i + N_w, I);
            j_min = max(j - N_w, 1);
            j_max = min(j + N_w, J);

            w_sum = 0;
            num   = 0;

            for jj = j_min:j_max
                for ii = i_min:i_max

                    pj2 = jj;
                    pi2 = ii;

                    patch_j = input_fillupBoundary(pi2:pi2+2*N_p, pj2:pj2+2*N_p);

                    d2 = sum((patch_j(:) - patch_i(:)).^2);
                    w_ij = exp(- d2 / denom);

                    w_ij = w_ij * (1 - abs(i-ii)/(N_w+1)) * (1 - abs(j-jj)/(N_w+1));

                    val = input(ii,jj);
                    num   = num   + w_ij * val;
                    w_sum = w_sum + w_ij;
                end
            end

            % normalize
            if w_sum > 0
                output(i,j) = num / w_sum;
            else
                output(i,j) = input(i,j);
            end
        end
    end
end