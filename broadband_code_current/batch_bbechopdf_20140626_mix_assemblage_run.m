% 2014 06 26  Run mixed assemblage model for comparison

mix_r = [5,10,30];
ns = 5e4;
gate_len = 0.5;
save_opt = 1;
tx = 3;
sdir = '/mnt/storage/ECHO_STAT/20140626_bbechopdf_mixed_assemblage';
if ~exist(sdir,'dir')
    mkdir(sdir);
end

Nweak = 1000;
N = [10,50,300];
for iR=1:length(mix_r)
    for iN=1:length(N)
        % Rayleigh
        bbechopdf([Nweak,N(iN)],[1,mix_r(iR)],ns,gate_len,tx,sdir,'rayleigh',[],0);
    end
end