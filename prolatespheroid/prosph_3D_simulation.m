function [data,varargout] = prosph_3D_simulation(e_ac,smpln,rough,varargin)
% Physics based simulation representing a randomly rough, randomly oriented
% prolate spheroid of eccentricity e_ac randomly located in a cylindrical
% aperture beampattern of ka = 2*pi.
%
% NOTE: The aspect ratio 'e_ac' here is the ratio of the minor axis (a)
% to the major axis (c). This is the inverse of the aspect ratio
% 'xi' defined in the thesis.

% Wu-Jung Lee 2012/10/18  Extend Saurav's 2D simulation to 3D
%             2014/11/23  Found out that theta_sph is actually
%                         pi/2-theta_wanted, have to correct this
%                         at some point

% Setting constants
r = 0.1;
c = ((r^3)/(e_ac^2))^(1/3); 
e_ba = 1;

if ~isempty(varargin)
    rayl_r = varargin{1};
else
    rayl_r = 1;
end

[theta_sph,phi_sph] = rand_rect_angle(smpln);

% Generating random variables
%theta_sph = rand(1, 5e6)*pi/2;
if rough==1
    roughness = raylrnd(ones(1,smpln)*rayl_r/sqrt(2));
else
    roughness = ones(1,smpln);
end

% Calculating echo amplitude
rsa = (((c^2)*(e_ba^2))/4)*((sin(atan(e_ac./tan(theta_sph)))./cos(theta_sph)).^2);

% Multiplying factors
data = rsa.*roughness;

% Binning data
[varargout{1}, varargout{2}] = logbinner(data, 300, 0);
%[pdf_x, pdf_y] = findEchoDist(data, 'default', 0,300);
