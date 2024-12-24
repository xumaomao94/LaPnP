function [C_match, nae, S] = Ctest(C_est,C,varargin)
%% match 
if ~all(C_est == 0)
    R = size(C,2);
    
    C_match = C_est; % K \times R
    for r = 1:R
        C_est(:,r) = C_est(:,r)/max(C_est(:,r));
        C(:,r) = C(:,r)/max(C(:,r));
    end
    
    % ind_match = 1:R;
    ind = [];
    for r = 1:R
        dist = sum( (C - kron(C_est(:,r),ones(1,R)) ).^2 ,1 );
        % dist = sum( (C_est - kron(C(:,r),ones(1,R)) ).^2 ,1);
        dist(ind) = inf;
        [~,ind(r)] = min(dist); % ind(r): C_est(:,r) - C(:,ind(r))
    end
    
    %$ compare
    C_match(:,ind) = C_match;
    % C_match = C_match(:,ind);
    
    nae = 0;
    for r = 1:R
        nae = nae+norm( C_match(:,r)/norm( C_match(:,r),1 ) - C(:,r)/norm(C(:,r)),1 );
    end
    nae = nae/R;
    
    S = 0;
    if ~isempty(varargin)
        S = varargin{1};
        S(:,:,ind) = S;
    end
else
    C_match = 0;
    nae = 0;
    S = 0;
end

end