% 2014 04 27  Figure out where the code was left in August 2013

mix_r = 1;
ns = 5e4;
gate_len = 0.5;
save_opt = 1;
tx = 3;
sdir = '/mnt/storage/ECHO_STAT/20140427_bbechopdf';
if ~exist(sdir,'dir')
    mkdir(sdir);
end

N = [2100:100:2400,2600:100:2900];
for iN=1:length(N)
    % Rayleigh
    bbechopdf(N(iN),mix_r,ns,gate_len,tx,sdir,'rayleigh',[],0);
end

N = [2100:100:2400,2600:100:2900];
for iN=1:length(N)
    % Fish, empirical length distribution
    bbechopdf(N(iN),mix_r,ns,gate_len,tx,sdir,'fish_defaultLenDistr_defaultAngleDistr',...
              [],1,[],[]);
end

%{
N = [10:10:90,100:100:900,1000,1500,2000,1100:100:1400,1600:100:1900,2500,3000,2100:100:2400,2600:100:2900];
for iN=1:length(N)
    % Prolate spheroid
    bbechopdf(N(iN),mix_r,ns,gate_len,tx,sdir,'prosph_0.2',...
              [],2,0.2);
end
%}
