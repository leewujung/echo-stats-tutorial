% 2014 04 27  Figure out where the code was left in August 2013
% 2014 05 26  Run more fish scatterer models for comparison

mix_r = 1;
ns = 5e4;
gate_len = 0.5;
save_opt = 1;
tx = 3;
sdir = '/mnt/storage/ECHO_STAT/20140626_bbechopdf_nobp';
if ~exist(sdir,'dir')
    mkdir(sdir);
end


N = [10,50,300,500,1000,2000];
for iN=1:length(N)
    % Rayleigh
    bbechopdf_nobp(N(iN),mix_r,ns,gate_len,tx,sdir,'rayleigh_nobp',[],0);
end