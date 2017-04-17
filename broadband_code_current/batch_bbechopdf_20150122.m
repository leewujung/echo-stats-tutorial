% 2014 04 27  Figure out where the code was left in August 2013
% 2014 05 26  Run more fish scatterer models for comparison
% 2014 09 04  Run bbechopdf with constraints on the polar angle
% 2015 01 22  Run fish scatterer models with +5deg swimbladder orientation

mix_r = 1;
ns = 5e4;
gate_len = 0.5;
save_opt = 1;
tx = 3;
sdir = '/mnt/storage/ECHO_STAT/20150122_fish_tiltedSwimbladder_bbechopdf';
if ~exist(sdir,'dir')
    mkdir(sdir);
end

N = [10:10:100,200:100:2000,2500:500:4000,2100:100:2400, 2600:100:2900, 3100:100:3400, 3600:100:3900];

for iN=1:length(N)
    % Fish, empirical length and angle of orientation distribution
    bbechopdf_20141124(N(iN),mix_r,ns,gate_len,tx,sdir,'fish_defaultLenDistr_defaultAngleDistr_tiltSwimbladder',...
              [],1,[],[-8,10]);
end
