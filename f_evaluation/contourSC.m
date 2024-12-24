function contourSC(mask,varargin)

    if iscell(varargin{1})
        varargin = varargin{1}; % allow cell input
    end
    SCcount = 0;
    SCindex = [];
    for j = 1:length(varargin)
        if ~isstring(varargin{j}) && ~iscell(varargin{j})
            SCcount = SCcount+1;
            SCindex = [SCindex,j];
        end
    end

    xmin = inf;
    xmax = -inf;
    for i = 1:SCcount %1:1 % 
        ind = SCindex(i);
        S = varargin{ind};
        valid_index = find(S>0);
        s = S(valid_index);
        xmin = min( xmin, min(s) );
        xmax = max( xmax, max(s) );
    end

    for i = 1:size(varargin{1},3)
        tiledlayout(1, SCcount, 'Padding', 'none', 'TileSpacing', 'compact'); 

        for j = 1:length(varargin)
            if isstring(varargin{j}) || iscell(varargin{j})
                title(varargin{j});
            else
                S = varargin{j};
                zero_index = find(mask==1);%(S<=1e-12);
                S = 10*log10(max(S,1e-23));
                S(zero_index) = NaN;

                nexttile
                contourf(kron(S,ones(2)), 200, 'linecolor', 'None');

                set(gca,'XTick',[],'YTick',[])
                colormap jet;
                clim([max(-140,10*log10(xmin(1))),min(10,10*log10(xmax(1)))+2])
                clim([-113,-70])
            end
        
        end
    end
    c = colorbar;
    c.Layout.Tile = 'east';
end