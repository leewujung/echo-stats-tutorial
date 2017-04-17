addpath /mnt/storage/prolatespheroid/
addpath /mnt/storage/modeling_code_current/

gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121129_midhann_width';

fname = 'midhann_narrow';
tx_opt = 1;  % square chirp
taper_opt = 7;  % mid Hann narrow
N = [50];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gate_len,tx_opt,...
                            taper_opt,save_opt,save_dir,fname);

end

fname = 'midhann_wide';
tx_opt = 1;  % square chirp
taper_opt = 9;  % mid Hann wide
N = [50];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gate_len,tx_opt,...
                            taper_opt,save_opt,save_dir,fname);

end

% No beampattern
gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121129_with_without_bp_bbmodel';
fname = 'tx_nosys';
tx_opt = 2;  % ideal tx, no system response
taper_opt = 0;
N = [50];
for iN=1:length(N)
    model_n_scat_no_bp(N(iN),sampleN,gate_len,tx_opt,...
                       taper_opt,save_opt,save_dir,fname);
end

