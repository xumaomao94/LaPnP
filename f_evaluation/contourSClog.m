function contourSClog(mask,varargin)

% e.g.,
% contourSC(X_true(:,:,48),"true map",XX{1}(:,:,48),"pnp",XX{2}(:,:,48),"interpolation",...
% XX{3}(:,:,48),"interMC",XX{4}(:,:,48),"Kriging",XX{5}(:,:,48),"Dowjons",XX{6}(:,:,48),"LL1")
    SCcount = 0;
    SCindex = [];
    for j = 1:length(varargin)
        if ~isstring(varargin{j})
            SCcount = SCcount+1;
            SCindex = [SCindex,j];
        end
    end

    xmin = inf;
    xmax = -inf;
    for i = 1:SCcount %1:1 % 
        ind = SCindex(i);
        S = varargin{ind};
        valid_index = find(S ~= -inf & S<0);
        s = S(valid_index);
        xmin = min( xmin, min(s) );
        xmax = max( xmax, max(s) );
    end
    % Linelevel = linspace(10*log10(xmin),10*log10(xmax),500);


    
    for i = 1:size(varargin{1},3)
        % figure(i)
        tiledlayout(SCcount, 1, 'Padding', 'none', 'TileSpacing', 'compact'); 
% tiledlayout(2, 4, 'Padding', 'none', 'TileSpacing', 'compact'); 
        % tiledlayout(2, 2, 'Padding', 'none', 'TileSpacing', 'compact'); 

        for j = 1:length(varargin)
            if isstring(varargin{j})
                title(varargin{j});
            else
                S = varargin{j};
                zero_index = find(mask==1);%(S<=1e-12);
                % S = 10*log10(max(S,1e-23));
                S(zero_index) = NaN;
                nexttile
                heatmap(S);
                % set(gca,'XTick',[],'YTick',[])
                colormap jet;
                clim([max(-110,xmin(1)),min(-40,xmax(1))])
                % clim([max(-90,10*log10(xmin(1))),min(10,10*log10(xmax(1)))])
                % clim([-105,-50])
                % clim([max(-140,10*log10(xmin(1))),-59])
                % clim([max(-140,10*log10(xmin)),-20])
                % clim([-130,-75])
                % clim([max(-13,10*log10(xmin)),min(10,10*log10(xmax))])
            end
        
        end
    end
    % c = colorbar;
    % c.Limits = [10*log10(xmin),10*log10(xmax)];
end