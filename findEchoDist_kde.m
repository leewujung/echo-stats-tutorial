% 2012 09 25  Use kernel density estimation to find echo pdf
%             Use Botev's kde method (Botev et al., 2010)

function [dens,x,bw] = findEchoDist_kde(s,npt)
% INPUT
%   s      samples
%   npt    number of points for echo pdf
% OUTPUT
%   dens   estimated density
%   x      x-axis of estimated density
%   bw     bandwidth of the kernel

% transform data into log10 space
[bw,dens_log,x_log]=kde(log10(s),npt,0);
%[dens_log,x_log,bw]=ksdensity(log10(s));
%x_log = x_log';
x = 10.^x_log.';   % transform back, adjust dimension as well
slope = 1./x;      % scale the estimated density back according
                   % to corresponding bandwidth on linear scale
scale = trapz(x,dens_log.*slope);  % rescale in the linear domain
dens = dens_log.*slope./scale;
