% 2013 07 27  beampattern response in the frequency domain

function bp = bpf_2way_fcn(theta,fvec,a)
% function [f,bp] = bpir_2way_fcn(theta,fvec,a)
% INPUT
%   theta   angles in [rad]
%   fvec    frequency vector [Hz]
%   a       radius of circular transducer [m]

if size(theta,1)~=1
    theta = theta.';
end
if size(fvec,2)~=1
    fvec = fvec.';
end

c = 1500;  % sound speed [m/s]
const = a/c*sin(theta);

bp = (2*besselj(1,fvec*const*2*pi)./(fvec*const*2*pi)).^2;  % two-way bp
bp(1,:) = 1;  % freq=0 case
bp(:,1) = 1;  % theta=0 case


