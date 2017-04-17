% 2011 12 08  Check again on the ir of beampattern response
%             because the previous one had a constant factor error in it
%             need to make sure everything is correct

function [t_conv,bp_t_conv] = bpir_2way_fcn(theta,fmax,a)
% function [t_conv,bp_t_conv] = bpir_2way_fcn(theta,fmax,a)
% INPUT
%   theta   angles in [rad]
%   fmax    max frequency, determines temporal resolution
%   a       radius of circular transducer [m]

%a = 0.054;  % [m]
c = 1500;  % sound speed [m/s]
%df = 1;
%fmax = 1e7;
%theta = 10/180*pi;  % polar angle from transducer axis
const = a/c*sin(theta);

%% Theoretical inverse fourier transform results
dt = 1/(2*fmax);
%dt = 1/(2*max(f));
t = -const:dt:const;
bp_t = 2/abs(const)*sqrt(2/pi)*sqrt(1-(t/const).^2); % 1-way bp, need
                                                     % convolution
bp_t_conv = conv(bp_t,bp_t)*dt/sqrt(2*pi);  % conv to get 2-way bp
             % sqrt(2*pi) is for the conversion between fft conventions

t_conv = ((1:length(bp_t_conv))-floor((length(bp_t_conv)+1)/2)-1)*dt;

%% Numerical ifft results
%f = 0:df:fmax;
%bp = (2*besselj(1,const*2*pi*f)./(const*2*pi*f)).^2;  % two-way bp
%bp(1) = 1;

%bp_full = [bp,fliplr(conj(bp(2:end)))];
%bp_ifft = length(bp_full)*ifftshift(ifft(bp_full))*df;
%bp_ifft = bp_ifft*sqrt(2*pi);  % sqrt(2*pi) is for the conversion between
%                               % differet forms of fourier transform
%t_ifft = ((1:length(bp_ifft))-length(bp)-1)*dt;
%bpt_half_len = floor((length(bp_full)+1)/2);

%% Plot
%figure;
%plot(t_ifft,abs(bp_ifft));
%hold on
%plot(t_conv,bp_t_conv,'r--');
%xlim([-10 10])

