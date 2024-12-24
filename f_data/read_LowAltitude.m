function S = read_LowAltitude(ssbIndex,GPSrange,M,N)
% GPSrange: [Longitude, Latitude] at
%           downmost(south)
%           leftmost(west)
%           rightmost(east)
%           top(north)

% map size: I * J, only support for a rectangular area


filename = ['4.9GHz PCI-131 at 98m relative height SSB-' num2str(ssbIndex) ' RSRP.mat'];
load(filename);

L = norm( GPSrange(3,:) - GPSrange(1,:) );
W = norm( GPSrange(2,:) - GPSrange(1,:) );

L_sep = L/M;
W_sep = W/N;

cos_theta = (GPSrange(3,1) - GPSrange(1,1))/L;
sin_theta = -(GPSrange(3,2) - GPSrange(1,2))/L;
Rmat_theta = [cos_theta,-sin_theta;sin_theta,cos_theta];

rsrp(:,[2,1]) = (Rmat_theta * rsrp(:,[2,1])')';
GPSrange = (Rmat_theta * GPSrange')';
% figure;
% A = 20 * ones(height(rsrp), 1);
% s1 = geoscatter(rsrp(:,1), rsrp(:,2), A, rsrp(:,3), "filled");
% c = colorbar;
% c.Label.String = "PBCH\_XSS\_RSRP";
% caxis([-120 -70]);



x = ceil( (rsrp(:,2) - GPSrange(1,1)) / L_sep );
x = min( max( 0,x ), M+1 );
y = ceil( (rsrp(:,1) - GPSrange(1,2)) / W_sep );
y = min( max( 0,y ), N+1);

S = zeros(M+2,N+2);
C = zeros(M+2,N+2);
for i = 1:height(rsrp)
    S(x(i)+1,y(i)+1) = rsrp(i,3) + S(x(i)+1,y(i)+1);
    C(x(i)+1,y(i)+1) = C(x(i)+1,y(i)+1) + 1;
end
S = S(2:end-1,2:end-1);
C = C(2:end-1,2:end-1);
mask = double(C ~= 0);

C(C==0) = 1;
S = S./C;
S = interp_random(S,mask,"NN");

% figure; contourf(kron(S',ones(3,3)), 200, 'linecolor', 'None'); colormap jet; set(gca,'XTick',[],'YTick',[]); clim([-120,-65]); c = colorbar
end