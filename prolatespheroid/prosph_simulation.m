function [pdf_x, pdf_y] = prosph_simulation(e_ac,smpln,rough)
% Physics based simulation representing a randomly rough, randomly oriented
% prolate spheroid of eccentricity e_ac randomly located in a cylindrical
% aperture beampattern of ka = 2*pi.

% WJL 2012 10 18  Add extra params for roughness and smpl number

% Setting constants
r = 0.1;
c = ((r^3)/(e_ac^2))^(1/3); 
e_ba = 1;

% Generating random variables
theta_sph = rand(1,smpln)*pi/2;
if rough==1
    roughness = raylrnd(ones(1,smpln)/sqrt(2));
else
    roughness = ones(1,smpln);
end

% Calculating echo amplitude
rsa = (((c^2)*(e_ba^2))/4)*((sin(atan(e_ac./tan(theta_sph)))./cos(theta_sph)).^2);

% Multiplying factors
data = rsa.*roughness;

% Binning data
[pdf_x, pdf_y] = logbinner(data, 300, 0);
%[pdf_x, pdf_y] = findEchoDist(data, 'default', 0,300);