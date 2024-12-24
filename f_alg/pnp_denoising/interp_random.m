function S = interp_random(S,Omega,varargin)
% fill the unobserved values, for 2d only

if isempty(varargin)
    interp_method = "NN";
else
    interp_method = varargin{1};
end

switch interp_method
    case "linear"
        % weight = [1 3 5 3 1; % linear interpolation
        %           3 5 7 5 3;
        %           5 7 9 7 5;
        %           3 5 7 5 3;
        %           1 3 5 3 1];
        weight = [1 3 1; % linear interpolation
                  3 5 3;
                  1 3 1];

        [I,J] = size(S);
        Ovec = find(Omega);
        
        % padding the boundary
        M = (size(weight,1)-1)/2;
        N = (size(weight,2)-1)/2;
    
        S_pad = zeros(I+2*M,J+2*N);
        S_pad(1+M:end-M,1+N:end-N) = S;
    
        O_count = zeros(I+2*M,J+2*N);
    
    
        for i = 1:length(Ovec)
            ind = Ovec(i);
            [x,y] = ind2sub([I,J],ind);
            
            S_pad(x:x+2*M,y:y+2*N) = S_pad(x:x+2*M,y:y+2*N) + S(x,y)*weight;
            O_count(x:x+2*M,y:y+2*N) = O_count(x:x+2*M,y:y+2*N) + weight;
        end
        
        S = (1-Omega) .* S_pad(1+M:end-M,1+N:end-N) ./ (O_count(1+M:end-M,1+N:end-N)+1e-9)...
            + Omega .* S;
        
    
        O_new = zeros(I,J);
        O_new( find( O_count(1+M:end-M,1+N:end-N) ) ) = 1;
    
        if sum( O_new(:) ) ~= I*J
            S = interp_random(S,O_new,"linear");
        end

    case "NN" % nearest neighbor
        % [I,J] = size(S);
        % Ovec_inv = find(1-Omega);
        % for i = 1:length(Ovec_inv)
        %     ind = Ovec_inv(i);
        %     [x,y] = ind2sub([I,J],ind);
        %     ind_match = find_neighbor(Omega,x,y);
        %     S(x,y) = mean(S(ind_match));
        % end
        weight = [1 1 1; % linear interpolation
                  1 1 1;
                  1 1 1];

        [I,J] = size(S);
        Ovec = find(Omega);
        
        % padding the boundary
        M = (size(weight,1)-1)/2;
        N = (size(weight,2)-1)/2;
    
        S_pad = zeros(I+2*M,J+2*N);
        S_pad(1+M:end-M,1+N:end-N) = S;
    
        O_count = zeros(I+2*M,J+2*N);
    
    
        for i = 1:length(Ovec)
            ind = Ovec(i);
            [x,y] = ind2sub([I,J],ind);
            
            S_pad(x:x+2*M,y:y+2*N) = S_pad(x:x+2*M,y:y+2*N) + S(x,y)*weight;
            O_count(x:x+2*M,y:y+2*N) = O_count(x:x+2*M,y:y+2*N) + weight;
        end
        
        S = (1-Omega) .* S_pad(1+M:end-M,1+N:end-N) ./ (O_count(1+M:end-M,1+N:end-N)+1e-9)...
            + Omega .* S;
        
    
        O_new = zeros(I,J);
        O_new( find( O_count(1+M:end-M,1+N:end-N) ) ) = 1;
    
        if sum( O_new(:) ) ~= I*J
            S = interp_random(S,O_new,"NN");
        end

    case "TPS"
        S = TPS(S,Omega);
        
end

end