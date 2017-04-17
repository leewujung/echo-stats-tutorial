% call model_n_scat_* in a batch


gate_len = 0.5;
sampleN = 1e4;
tx_opt = 3;  % actual transmit, with system response
taper_opt = 0;  % no tapering

save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121101_bb_smpl1e4_noise';

addpath /mnt/storage/prolatespheroid/
addpath /mnt/storage/modeling_code_current/

fname = 'tx_sys';
N = [450:50:1000,1500,2000];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gate_len,tx_opt,...
                            taper_opt,save_opt,save_dir,fname);
end
