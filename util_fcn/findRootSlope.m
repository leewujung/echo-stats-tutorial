function [root,slope] = findRootSlope(b,xs,t,dt_f,ka)

t_fine = t(xs-1):dt_f:t(xs+1); % theta
bt_fine = (2*besselj(1,ka*sin(t_fine))./(ka*sin(t_fine))).^2; %b(theta)
bt_fine = real(bt_fine);
[dum,pos] = min(abs(b-bt_fine));
root = t_fine(pos);
slope = (bt_fine(pos+1)-bt_fine(pos-1))/(2*dt_f);
% slope = (bt_fine(pos)-bt_fine(pos-1))/dt_f;


