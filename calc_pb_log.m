% 2010 12 07
function [b,pb] = calc_pb_log(ka,b_start_log,b_end_log,b_num)
% function [b,pb] = calc_pb(ka,del_b)
%
% Calculate the beampattern pdf Pb(b)
%
% ka       ka value for piston transducer
% del_b    resolution of the resultant pb
%          this value is not related to the analytical part in the beginning
%          of this function
%
% WJL 2011/03/26
% 2016 07 18  small edit for locs empty case
%             remove wait bar
% 2016 10 18  change to log spacing

% ka=2*pi;

%% Calculate beampattern pdf Pb(b)
% analytical solution ===========================================
dtheta = 1e-6*pi;
theta = 0:dtheta:pi/2;
theta(1) = 1e-16;
btheta = (2*besselj(1,ka*sin(theta))./(ka*sin(theta))).^2;
% btheta(1) = 1;
diffb = diff(btheta);

[pks,locs] = findpeaks(btheta);
if isempty(locs)
    include = 1:length(btheta);
else
    include = find(btheta>btheta(locs(1))); % discard multi-valued part
end

theta_inc = theta(include);
pb_analy = sin(theta_inc).^2./...
           (4*cos(theta_inc).*besselj(2,ka*sin(theta_inc)).*sqrt(btheta(include)));
btheta_analy = btheta(include);


% numerical solution ============================================
dtheta = 1e-5*pi;
theta = 0:dtheta:pi/2;
btheta = (2*besselj(1,ka*sin(theta))./(ka*sin(theta))).^2;
btheta(1) = 1;

% del_b = 0.00115;
% del_b = 1e-6;
dtheta_fine = 1e-6;
% b = 0:del_b:1;  % original case before 2016/10/18
b = logspace(b_start_log,b_end_log,b_num);
% b = logspace(-3,0,1e4);
% b(1) = [];
% b = [logspace(-4,-3,1e2),b];

pb = zeros(1,length(b)-1);
%h = waitbar(0,'Please wait...');
for iB=1:(length(b)-1)
    %waitbar(iB/(length(b)-1));
    x = b(iB)- btheta;
    xs = (find(diff(sign(x))~=0)+1)';
    xs3 = [xs-1 xs xs+1];
    [dum,xs3pos] = min(abs(x(xs3)),[],2);
    x_pos = xs-2+xs3pos; % zero-crossing position

    root = zeros(1,length(x_pos));
    slope = zeros(1,length(x_pos));
    for iR=1:length(x_pos)
        [root(iR),slope(iR)] = findRootSlope(b(iB),x_pos(iR),theta,dtheta_fine,ka);
    end
    pb(iB) = sum(sort(sin(root)./abs(slope)));
end
%close(h);
pb(length(pb)+1) = pb_analy(1); % take the last value (the maximum at theta=0)
                                % using analytical soln



