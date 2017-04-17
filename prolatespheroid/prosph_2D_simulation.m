function [pdf_x, pdf_y, data] = prosph_2D_simulation(e_ac,smpln,rough,varargin)
% Physics based simulation representing a randomly rough, randomly oriented
% prolate spheroid of eccentricity e_ac randomly located in a cylindrical
% aperture beampattern of ka = 2*pi.

% varargin{1}  noise rms
% varargin{2}  range of angles: in terms of the fraction of pi/2

% WJL 2012 10 18  Add extra params for roughness and smpl number
% WJL 2012 12 02  Flexible to add noise
% WJL 2013 01 14  Allow arbitrary range of angles of orientation

if ~isempty(varargin)
    if length(varargin)==1
        fnrms = varargin{1};
        a_range = 1;
    else
        fnrms = varargin{1};
        a_range = varargin{2};
    end
else
    fnrms = [];
    a_range = 1;
end

% Setting constants
r = 0.1;
c = ((r^3)/(e_ac^2))^(1/3); 
e_ba = 1;

% Generating random variables
theta_sph = rand(1,smpln)*pi/2*a_range;
if rough==1
    roughness = raylrnd(ones(1,smpln)/sqrt(2));
else
    roughness = ones(1,smpln);
end

% Calculating echo amplitude
rsa = (((c^2)*(e_ba^2))/4)*((sin(atan(e_ac./tan(theta_sph)))./cos(theta_sph)).^2);

% Multiplying factors
data = rsa.*roughness;

drms = sqrt(mean(data.^2));
if ~isempty(fnrms)
    n = normrnd(0,1/sqrt(2)/fnrms*drms,size(data));
else
    n = zeros(size(data));
end
%data = abs(data + n);
data = data + n;
data(data<0) = [];

% Binning data
[pdf_x, pdf_y] = logbinner(data, 300, 0);
%[pdf_x, pdf_y] = findEchoDist(data, 'default', 0,300);