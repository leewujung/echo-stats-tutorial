

N = [90,100:100:1000,1200,1500,2000];
%N = 10;
mix_r = [1];
sampleN = 5e4;
gate_len = 0.5;
save_opt = 1;
taper_opt = 0;  % taper transmit signal
tx_opt = 3;
fname = 'tx_sys';
save_dir0 = '/mnt/storage/ECHO_STAT/20130801_bb_pdfmodel_rayleigh';
indiv0 = 0;
save_dir1 = '/mnt/storage/ECHO_STAT/20130801_bb_pdfmodel_fish';
indiv1 = 1;

for iN=1:length(N)
    model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,sampleN,gate_len,tx_opt, ...
                                             taper_opt,save_opt,save_dir0,fname,indiv0);
    model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,sampleN,gate_len,tx_opt, ...
                                             taper_opt,save_opt,save_dir1,fname,indiv1);
end
