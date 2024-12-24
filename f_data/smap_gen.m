function [X,Sc,Ctrue] = smap_gen(map_size, R, shadow_sigma, d_corr, varargin)

%%  parse parameters
p = inputParser;
addOptional(p,'basis',"sinc",@isstring); % "gaussian": gaussian kernel
addOptional(p,'num_peaks_per_psd',3,@isnumeric);
addOptional(p,'strictly_separable',false,@islogical); % true: peaks are generated independently for each emitter
addOptional(p,'seed',sum(100*clock),@isnumeric); % set the same so all methods perform on the same map
parse(p,varargin{:});
par = p.Results;

basis = par.basis;
num_peaks_per_psd = par.num_peaks_per_psd;
strictly_separable = par.strictly_separable;
seed = par.seed;
rng(seed)

K = map_size(3);

%% generate psd of the emitters
indK = (1:K)';
if strcmp(basis, "gaussian")
    psd = @(f0,sigma) exp(-(indK-f0).^2 /(2*sigma^2)); 
elseif strcmp(basis, "sinc")
    psd = @(f0,psd_span) sinc((indK-f0)/psd_span).^2.*( abs((indK-f0)/psd_span)<=1);
end
    

if strictly_separable
    ind_psd = R*3+2:2:K-1; % the first pulses are defined as 4,6,..
else
    ind_psd = 2:K-1; 
end

Ctrue = zeros(K,R);
if K ~= 1
    if ~strictly_separable
        ind_psd = 2:2:K-2;
        for rr = 1:R
            psd_peaks = ind_psd(   randperm(length(ind_psd), num_peaks_per_psd)   );
            am = 0.5 + 1.5*rand(num_peaks_per_psd,1);
    
            cr = 0;
            for q = 1:num_peaks_per_psd
                cr = cr + am(q)*psd(   psd_peaks(q),  4+4*rand());
            end
            Ctrue(:,rr) = cr/norm(cr);
        end
    else
        for rr = 1:R
            psd_peaks = ind_psd(   randperm(length(ind_psd), num_peaks_per_psd)   );
            am = 0.5 + 1.5*rand(num_peaks_per_psd,1);
    
            % First peak, ensure separability if needed
            cr = am(1) * psd(   1+(rr-1)*3,   1.5 + 1.5*rand()); % 
    
            for q = 2:num_peaks_per_psd
                cr = cr + am(q)*psd(   psd_peaks(q),   2+2.5*rand());
            end
            Ctrue(:,rr) = cr/norm(cr);
        end
    end
else
    for rr = 1:R
        Ctrue(:,rr) = 0.5 + 1.5*rand(1,1);
    end
end

%% generate spatial map
loss_f = @(x,d,alpha) min(1,(x/d).^(-alpha));
d0 = 2.5;

[Xmesh_grid, Ymesh_grid] = meshgrid(0:map_size(1)-1, 0:map_size(2)-1);
Xgrid = Xmesh_grid + 1i*Ymesh_grid;

Sc = zeros(map_size([1,2]));
for rr = 1:R
    location = map_size(1)*rand() + 1i*map_size(1)*rand();
    peak = [real(location), imag(location)];
    
    loss_mat = abs(Xgrid - location);
    alpha = 2+0.5*rand();
    p = exp(-1/d_corr);
    shadow = Shadowing(Xgrid,shadow_sigma,p);
    shadow_linear = 10.^(shadow/10);
    
    Sc_r = loss_f(loss_mat,d0,alpha).*shadow_linear;
    Sc_r = Sc_r/norm(Sc_r,'fro');
    Sc(:,:,rr) = Sc_r;
end

Sc_vec = reshape(Sc,[],R);
X = Sc_vec * Ctrue.';
X = reshape(X,map_size);

end