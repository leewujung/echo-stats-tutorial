

N = [10,20,50,100];
%N = [10:10:90,100:100:2000];
%N = 10;
mix_r = 1;
ns = 5e4;
gate_len = 0.5;
save_opt = 1;
tx = 3;
sdir = '/mnt/storage/ECHO_STAT/20130808_rayleigh_old_new_test';
if ~exist(sdir,'dir')
    mkdir(sdir);
end

for iN=1:length(N)

% bbechopdf(N,mix_r,ns,gate_len,tx,sdir,fname,...
%           bpa,indiv,[len_bin,len_dist],[angle_mean,angle_std],nb_freq);

% Rayleigh
bbechopdf(N(iN),mix_r,ns,gate_len,tx,sdir,'rayleigh',[],0);
model_n_scat_freqdep_bp_mixed_freqdomain(N(iN),mix_r,ns,gate_len,3,0,...
                                         1,sdir,'rayleigh_old_routine');

end