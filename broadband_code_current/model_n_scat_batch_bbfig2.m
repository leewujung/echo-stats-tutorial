% 2012 11 08  Calculate echo pdf model


addpath /mnt/storage/prolatespheroid/
addpath /mnt/storage/modeling_code_current/

% monotype aggregation, no system response
gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121128_tx_nosys_bbmodel_2ndset';
fname = 'tx_nosys';
tx_opt = 2;  % ideal tx
taper_opt = 0;  % no additional taper
%N = [10,50,300,500:100:1000,20:10:40,60:10:90,100,200,400];
N = [50];
r = 1;
for iN=1:length(N)
    model_n_scat_freqdep_bp_mixed(N(iN),r,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,0);
end

% monotype aggregation, with system response
gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121128_pdfmodel_sample_only_2ndset';
fname = 'tx_sys';
tx_opt = 3;  % actual tx
taper_opt = 0;  % no additional taper
%N = [10,50,300,500:100:1000,20:10:40,60:10:90,100,200,400];
N = [10,300,50];
r = 1;
for iN=1:length(N)
    model_n_scat_freqdep_bp_mixed(N(iN),r,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,0);
end
