function ind_match = find_neighbor(O,x,y,varargin)
    % find the nearest neighbor of (x,y) in O

    if isempty(varargin)
        d = 1;
    else
        d = varargin{1};
    end

    [I,J] = size(O);
    x_min = x - d;
    x_max = x + d;
    y_min = y - d;
    y_max = y + d;

    search_range = [4*d,2];
    for i = 1:d
        search_range((i-1)*4+1 : i*4, :) = [
                        [x_min+i, y+i];
                        [x+i,y_max-i];
                        [x_max-i,y-i];
                        [x-i,y_min+i]];
    end
    for i = 1:4*d
        if search_range(4*d-i+1,1) <=0 || search_range(4*d-i+1,2) <= 0 ||...
            search_range(4*d-i+1,1) >I || search_range(4*d-i+1,2) > J
            search_range(4*d-i+1,:) = [];
        end
    end
    
    search_ind = sub2ind([I,J],search_range(:,1),search_range(:,2));

    O_temp = zeros([I,J]);
    O_temp(search_ind) = O(search_ind);
    
    ind_match = find( O_temp );
    if isempty(ind_match)
        ind_match = find_neighbor(O,x,y,d+1);
    end

end