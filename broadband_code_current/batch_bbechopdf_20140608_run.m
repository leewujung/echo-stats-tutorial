% 2014 04 27  Figure out where the code was left in August 2013
% 2014 05 26  Run more fish scatterer models for comparison

mix_r = 1;
ns = 5e4;
gate_len = 0.5;
save_opt = 1;
tx = 3;
sdir = '/mnt/storage/ECHO_STAT/20140608_bbechopdf';
if ~exist(sdir,'dir')
    mkdir(sdir);
end


%N = [10:10:100,200:100:2000,2500,3000];
%N = [10,20,50,100,200,500,1000];
%N = [2100:100:2400,2600:100:2900];
%{
for iN=1:length(N)
    % Fish, empirical length distribution
    bbechopdf(N(iN),mix_r,ns,gate_len,tx,sdir,'fish_defaultLenDistr_mean5_std_10',...
              [],1,[],[5,10]);
end
%}

N = [100];
for iN=1:length(N)
    % Fish, empirical length distribution
    bbechopdf(N(iN),mix_r,ns,gate_len,tx,sdir,'fish_defaultLenDistr_mean13_std_10',...
              [],1,[],[13,10]);
    bbechopdf(N(iN),mix_r,ns,gate_len,tx,sdir,'fish_defaultLenDistr_mean-13_std_10',...
              [],1,[],[-13,10]);
end


%{
N = [10:10:90,100:100:900,1000,1500,2000,1100:100:1400,1600:100:1900,2500,3000,2100:100:2400,2600:100:2900];
for iN=1:length(N)
    % Prolate spheroid
    bbechopdf(N(iN),mix_r,ns,gate_len,tx,sdir,'prosph_0.2',...
              [],2,0.2);
end
%}
