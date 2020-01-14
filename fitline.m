function [ xy,r, alpha ] = fitline( xyin )
[alpha r] = get_alphaR(xyin(:,1), xyin(:,2));
theta_line=[0:0.001: 2*pi]';
rho = r./(cos(theta_line-alpha));
d=rho.*cos(theta_line-alpha)-r;
xy=[cos(theta_line).*rho,sin(theta_line).*rho];


% Old End Finding method:
% idxx=find(xy(:,1)>min(xyin(:,1)) & xy(:,1)<max(xyin(:,1)));
% idxy=find(xy(:,2)>min(xyin(:,2)) & xy(:,2)<max(xyin(:,2)));
% idx=intersect(idxx,idxy);
% xy=[xy(idx,1),xy(idx,2)];


% Method to find ends of line:
distancesstart=sqrt((xyin(1,1)-xy(:,1)).^2+(xyin(1,2)-xy(:,2)).^2);
idxstart=find(distancesstart==min(distancesstart));
distancesend=sqrt((xyin(end,1)-xy(:,1)).^2+(xyin(end,2)-xy(:,2)).^2);
idxend=find(distancesend==min(distancesend));
xy=[xy(idxstart,1),xy(idxstart,2);
xy(idxend,1), xy(idxend,2)];
end
