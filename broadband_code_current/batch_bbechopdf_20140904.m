% 2014 04 27  Figure out where the code was left in August 2013
% 2014 05 26  Run more fish scatterer models for comparison
% 2014 09 04  Run bbechopdf with constraints on the polar angle

mix_r = 1;
ns = 5e4;
gate_len = 0.5;
save_opt = 1;
tx = 3;
sdir = '/mnt/storage/ECHO_STAT/20140904_bbechopdf_restricted_polar_angle';
if ~exist(sdir,'dir')
    mkdir(sdir);
end


N = [10,20,50,100];

%{
for iN=1:length(N)
    % Rayleigh scatterer
    bbechopdf(N(iN),mix_r,ns,gate_len,tx,sdir,'bpa_30deg_rayleigh',...
              30,0);
end
%}

for iN=1:length(N)
    % Rayleigh scatterer
    bbechopdf(N(iN),mix_r,ns,gate_len,tx,sdir,'rayleigh',...
              [],0);
end

for iN=1:length(N)
    % Fish, empirical length and angle of orientation distribution
    bbechopdf(N(iN),mix_r,ns,gate_len,tx,sdir,'bpa_30deg_fish_defaultLenDistr_defaultAngleDistr',...
              30,1,[],[-13,10]);
end

for iN=1:length(N)
    % Prolate spheroid
    bbechopdf(N(iN),mix_r,ns,gate_len,tx,sdir,'bpa_30deg_prosph_0.5',...
              30,2,0.5);
end

