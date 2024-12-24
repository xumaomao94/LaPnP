function [output,W] = DSG_NLM(input, sqrNoise, N_w, N_p, varargin)

% N_w: search window size = 2 N_w + 1
% N_p: patch size = 2 N_p + 1

    Npower = sqrNoise^2;
    patch_size = 2*N_p + 1;
    window_size = 2*N_w + 1;

    [I,J] = size(input); % only support 2D currently
    N = I*J;

    W = zeros(N,N);
    % W = sparse(W);

    input_fillupBoundary = zeros(I+2*N_p,J+2*N_p);
    input_fillupBoundary(N_p+1:N_p+I, N_p+1:N_p+J) = input;
    input_fillupBoundary(N_p+1:N_p+I, 1:N_p) = kron(input(:,1),ones(1,N_p));
    input_fillupBoundary(N_p+1:N_p+I, N_p+J+1:2*N_p+J) = kron(input(:,end),ones(1,N_p));
    input_fillupBoundary(1:N_p,:) = ...
        kron(input_fillupBoundary(N_p+1,:),ones(N_p,1));
    input_fillupBoundary(N_p+I+1:2*N_p+I,:) = ...
        kron(input_fillupBoundary(N_p+I,:),ones(N_p,1));
    
    input_fillupBoundary = log10(1e-15 + input_fillupBoundary);
    %% Filter weights
    if isempty(varargin)
        for i = 1:N
            % [m,n]-th pixel
            n = ceil(i/I);
            m = i - (n-1)*I;
    
            % patch: m+Np-Np:m+Np+Np, n+Np-Np:n+Np+Np
            patch_i = input_fillupBoundary( m:m+2*N_p,n:n+2*N_p );
            % norm_i = sum( (patch_i(:)+1e-12).^2 );

            % search window: (m_min,n_min; m_max,n_max)
            search_Window = [max(m-N_w,1), max(n-N_w,1);
                             min(m+N_w,I), min(n+N_w,J)];
    
            for nn = n : search_Window(2,2)
                if nn == n
                    mm_start = m;
                else
                    mm_start = search_Window(1,1);
                end
                
                for mm = mm_start:search_Window(2,1)
                    j = (nn-1)*I + mm;
                    patch_j = input_fillupBoundary( mm:mm+2*N_p,nn:nn+2*N_p );
                    
                    % W(i,j) = exp( - sum((patch_j(:) - patch_i(:)).^2) / norm_i / 2  / Npower);
                    W(i,j) = exp( - sum((patch_j(:) - patch_i(:)).^2) / patch_size^2 / 2  / Npower);
                   

                    W(i,j) = W(i,j) * (1 - abs(m-mm)/(N_w+1)) * (1 - abs(n-nn)/(N_w+1));
                end
            end
        end
        
        W = W + triu(W).' - diag(diag(W));
    
        % Normalize
        w_sum_1 = sum(W,1); % sum all rows, 1*N
        w_sum_2 = sum(W,2); % sum all columns, N*1
    
        W = W ./ sqrt( w_sum_2*w_sum_1 );
        
        w_sum_1 = sum(W,1);
        alpha = 1/max(w_sum_1);
        W = alpha*W;
        w_sum_1 = sum(W,1);
        W = W + diag( 1-w_sum_1 );
    else % input W only, otherwise it would be wrong
        W = varargin{1};
    end
    
    %% denoise
    output = W * input(:);
    output = reshape(output,[I,J]);
end