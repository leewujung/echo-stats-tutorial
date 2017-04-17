
addpath /mnt/storage/prolatespheroid/
addpath /mnt/storage/modeling_code_current/

%{
% monotype aggregation
gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121108_pdfmodel_sample_only';
fname = 'tx_sys';
tx_opt = 3;  % actual tx
taper_opt = 0;  % no additional taper
%N = [10,50,300,500:100:1000,20:10:40,60:10:90,100,200,400];
N = [3200:200:4000];
r = 1;
for iN=1:length(N)
    model_n_scat_freqdep_bp_mixed(N(iN),r,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,0);
end
%}

% monotype aggregation
gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121116_pdfmodel_w_respenv_5e4';
fname = 'tx_sys';
tx_opt = 3;  % actual tx
taper_opt = 0;  % no additional taper
%N = [10,50,300,500:100:1000,20:10:40,60:10:90,100,200,400];
N = [1200:200:3000];
r = 1;
for iN=1:length(N)
    model_n_scat_freqdep_bp_mixed(N(iN),r,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,0);
end

%{
% monotype aggregation
gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121127_tx_nosys_bbmodel';
fname = 'tx_nosys';
tx_opt = 2;  % ideal tx
taper_opt = 0;  % no additional taper
N = [20:10:40,60:10:90,100,200,400:100:1000,1500,2000];
r = 1;
for iN=1:length(N)
    model_n_scat_freqdep_bp_mixed(N(iN),r,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,0);
end
%}