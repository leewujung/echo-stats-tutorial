function [data,varargout] = prosph_amp(e_ac,smpln,rough,varargin)

% INPUT
% e_ac        ratio of the minor axis to major axis
% smpln       desired number of samples
% rough       =1: rough, =0: smooth
% varagin{1}  ratio of strong to weak scatterer
% varagin{2}  [mean,std] of angle of orientation distribution
% varagin{3}  [len_bin,len_dist] length distribution

% Wu-Jung Lee 2012 10 18  Extend Saurav's 2D simulation to 3D
%             2014 09 04  modify prosph_3D_simulation to take angle
%                         of orientation distribution
%                         also add in length variation
%             2014 11 23  realized that the theta_sph generated
%                         using rand_rect_angle.m actually
%                         corresponds to pi/2-theta_sph. Therefore
%                         the results generated assuming all angles
%                         of orientation prior to this date is
%                         wrong and has to be calculated again
%             2014 11 24  modify the code so that theta_sph can use
%                         [0,2*pi] uniform distribution when
%                         [angle_mean,angle_std]=[NaN,NaN]

% Setting constants
%r = 0.1;
%c = ((r^3)/(e_ac^2))^(1/3); 
%e_ba = 1;

if nargin>3
    if ~isempty(varargin{1})
        rayl_r = varargin{1};
    else
        rayl_r = 1;
    end
end

if nargin>4
    if ~isempty(varargin)  % with specified angle of orientation
        angle_mean = varargin{2}(1);
        angle_std = varargin{2}(2);
        if ~isnan(angle_mean)
            disp('with assigned angle of orientation');
            theta_sph = normrnd_truncated(angle_mean,angle_std,2,smpln,[]); % truncated normal distribution
        else
            theta_sph = unifrnd(zeros(1,smpln),2*pi); % default to [0,2*pi] uniform distribution
        end
    else  % default to [0,2*pi] uniform distribution
        theta_sph = unifrnd(zeros(1,smpln),2*pi);
    end
else  % default to [0,2*pi] uniform distribution
    theta_sph = unifrnd(zeros(1,smpln),2*pi);
end

if nargin>5
    if ~isempty(varargin{3})
        len_bin = varargin{3}(:,1);
        len_dist = varargin{3}(:,2);
    else
        disp('Use default fish length distribution');
        fishL = load(['/mnt/storage/broadband_code_current/fish_info/fish_len_dist.mat']);
        len_bin = fishL.L_bin;
        len_dist = fishL.L_dist;
    end
else
    disp('Use default fish length distribution');
    fishL = load(['/mnt/storage/broadband_code_current/fish_info/fish_len_dist.mat']);
    len_bin = fishL.L_bin;
    len_dist = fishL.L_dist;
end

% Generate random length
[len,~] = discrete_rnd(len_bin,len_dist,smpln);
cc = len'/2;   % length of semi-major axis
b1 = cc*e_ac; % length of semi-minor axis

% Calculating echo amplitude
fss = cc/2.*sin(atan(b1./(cc.*tan(theta_sph)))).^2./cos(theta_sph).^2;

% Generate roughness
if rough==1
    roughness = raylrnd(ones(1,smpln)*rayl_r/sqrt(2));
else
    roughness = ones(1,smpln);
end

% Multiplying factors
data = fss.*roughness;

% Binning data
[varargout{1}, varargout{2}] = logbinner(data, 300, 0);
