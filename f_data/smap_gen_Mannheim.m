function [X,mask] = smap_gen_Mannheim(varargin)

%%  parse parameters
p = inputParser;

default_filename = strcat(pwd,'\DATA\Manheim\formatted_data.mat');
addOptional(p,'data_file',default_filename,@isstring);

parse(p,varargin{:});
par = p.Results;


%% Obtain real dataset
load(par.data_file);

T = permute(T,[2,1,3]);

mask = sum(T(:,:,1:9),3) == 0;

X = 10.^(T(:,:,1:9)./10);
for k = 1:9
    Xk = X(:,:,k);
    Xk(mask) = 0;
    X(:,:,k) = Xk;
end

mask = double(mask);
end