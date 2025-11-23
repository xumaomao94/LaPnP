function [X,Sc,Ctrue,mask] = smap_gen_RMSeer(K,varargin)

%%  parse parameters
p = inputParser;
addOptional(p,'basis',"sinc",@isstring); % "gaussian": gaussian kernel
addOptional(p,'num_peaks_per_psd',3,@isnumeric);
addOptional(p,'strictly_separable',false,@islogical); % true: peaks are generated independently for each emitter
addOptional(p,'seed',sum(100*clock),@isnumeric); % set the same so all methods perform on the same map
addOptional(p,'RTmethod',"IRT2",@isstring); % the RT method in RadioMapSeer
addOptional(p,'map_number',0,@isnumeric);
addOptional(p,'emitter_number',[0,1],@isnumeric); % which emitter to choose, set 0 to choose all
addOptional(p,'downsample_rate',0.5,@isnumeric);

parse(p,varargin{:});
par = p.Results;

basis = par.basis;
num_peaks_per_psd = par.num_peaks_per_psd;
strictly_separable = par.strictly_separable;
seed = par.seed;
% rng('default')
rng(seed)

%% Obtain RadioMapSeer Smap
building_folder = fullfile(pwd, 'DATA', 'RadioMapSeer', 'png', 'buildings_complete');
building_file = fullfile(building_folder, sprintf('%d.png', par.map_number));

mask = imread(building_file);
mask = double(mask)./255;
mask = imresize(mask, par.downsample_rate);
mask(mask>0.1) = 1;
mask(mask<=0.1) = 0;

Sc_folder = fullfile(pwd, 'DATA', 'RadioMapSeer', 'gain');
R = length(par.emitter_number);
Sc = zeros([256*par.downsample_rate,256*par.downsample_rate,R]);
for r = 1:R
    Sc_file = fullfile(Sc_folder, par.RTmethod, ...
        sprintf("%d_%d.png", par.map_number,par.emitter_number(r)));
    Sc_temp = double(imread(Sc_file))./255;
    Sc(:,:,r) = imresize(Sc_temp, par.downsample_rate);
end

PLmax = -47;
PLth = -127;
PLtrunc = -147;

Sc = Sc * 100 + (-147); % -147 - -47
Sc = 10.^(Sc./10);
for r = 1:R
    Sc(:,:,r) = Sc(:,:,r).*(1-mask);
end

S_size = size(Sc);
map_size = [S_size(1:2),K];
R = S_size(3);

%% generate psd of the emitters
indK = (1:K)';
if strcmp(basis, "gaussian")
    psd = @(f0,sigma) exp(-(indK-f0).^2 /(2*sigma^2)); 
elseif strcmp(basis, "sinc")
    psd = @(f0,psd_span) sinc((indK-f0)/psd_span).^2.*( abs((indK-f0)/psd_span)<=1);
end
    

if strictly_separable
    ind_psd = R*3+3:1:K; % the first pulses are defined as 4,6,..
else
    % ind_psd = R*3+4:2:K-2; % the first pulses are defined as 4,6,..
    ind_psd = 2:K-1; 
end

Ctrue = zeros(K,R);
if K ~= 1
    if ~strictly_separable
        for rr = 1:R
            psd_peaks = ind_psd(   randperm(length(ind_psd), num_peaks_per_psd)   );
            am = 0.5 + 1.5*rand(num_peaks_per_psd,1);
    
            cr = 0;
            for q = 1:num_peaks_per_psd
                cr = cr + am(q)*psd(   psd_peaks(q),   4+4*rand());
            end
            Ctrue(:,rr) = cr/norm(cr);
        end
    else
        for rr = 1:R
            psd_peaks = ind_psd(   randperm(length(ind_psd), num_peaks_per_psd)   );
            am = 0.5 + 1.5*rand(num_peaks_per_psd,1);
    
            cr = am(1) * psd(   1+(rr-1)*3,   1.5+1.5*rand()   );
            
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

Sc_vec = reshape(Sc,[],R);
X = Sc_vec * Ctrue.';
X = reshape(X,map_size);

end