function [X,S,Ctrue] = smap_gen_chen(K,ssbIndex,varargin)
    
%%  parse parameters
p = inputParser;
addOptional(p,'M',52,@isnumeric);
addOptional(p,'N',60,@isnumeric);
addOptional(p,'basis',"sinc",@isstring); % "gaussian": gaussian kernel
addOptional(p,'num_peaks_per_psd',3,@isnumeric);
addOptional(p,'strictly_separable',true,@islogical); % true: peaks are generated independently for each emitter
addOptional(p,'seed',sum(100*clock),@isnumeric); % set the same so all methods perform on the same map

parse(p,varargin{:});
par = p.Results;
M = par.M;
N = par.N;
R = length(ssbIndex);

basis = par.basis;
num_peaks_per_psd = par.num_peaks_per_psd;
strictly_separable = par.strictly_separable;
seed = par.seed;
% rng('default')
rng(seed)

%% load
load("vec.mat");
GPSrange = [GPSvec(60,[2,1]);
            GPSvec(1,[2,1]);
            GPSvec(3120,[2,1]);
            GPSvec(3061,[2,1])];

S = zeros(M,N,length(ssbIndex));
for r = 1:length(ssbIndex)
    S(:,:,r) = read_LowAltitude(ssbIndex(r),GPSrange,M,N);
end
S = permute(S,[2,1,3]);
S = 10.^(S./10);

%% generate psd of the emitters
indK = (1:K)';
if strcmp(basis, "gaussian")
    psd = @(f0,sigma) exp(-(indK-f0).^2 /(2*sigma^2)); 
elseif strcmp(basis, "sinc")
    psd = @(f0,psd_span) sinc((indK-f0)/psd_span).^2.*( abs((indK-f0)/psd_span)<=1);
end
    

ind_psd = 2:2:K-2;%R*3+4:2:K-2; % the first pulses are defined as 4,6,..

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


S_vec = reshape(S,[],R);
X = S_vec * Ctrue.';
X = reshape(X,[N,M,K]);


S = S(4:54,1:51,:);
X = X(4:54,1:51,:);

end