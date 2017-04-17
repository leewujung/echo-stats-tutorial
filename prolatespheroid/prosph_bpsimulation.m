function [pdf_x, pdf_y] = prosph_bpsimulation(e_ac)
% Physics based simulation representing a randomly rough, randomly oriented
% prolate spheroid of eccentricity e_ac randomly located in a cylindrical
% aperture beampattern of ka = 2*pi.

% Setting constants
r = 0.1;
c = ((r^3)/(e_ac^2))^(1/3); 
e_ba = 1;
ka = 2*pi;

% Generating random variables
% x_st = linspace(1e-5, pi/2, 1e5);
% cdf_st = 1 - cos(x_st);
% theta = randarb(x_st, cdf_st, 5e6);
theta = rand(1, 1e7)*pi/2;
theta_sph = rand(1, 1e7)*pi/2;
roughness = raylrnd(ones(1, 1e7)/sqrt(2));

% Calculating beampattern
nump = 2*besselj(1, ka*sin(theta));
denp = 2*pi*sin(theta);
bp = (nump./denp).^2;

% Calculating echo amplitude
rsa = (((c^2)*(e_ba^2))/4)*((sin(atan(e_ac./tan(theta_sph)))./cos(theta_sph)).^2);

% Multiplying factors
data = rsa.*roughness.*bp;

% Binning data
[pdf_x, pdf_y] = logbinner(data, 300, 10);
% [pdf_y, pdf_x] = ksdensity(data, logspace(log10(min(data)) + 10, log10(max(data)), 300), 'width', 2e-6);
[pdf_x, pdf_y] = pdf_normalizer(pdf_x, pdf_y);