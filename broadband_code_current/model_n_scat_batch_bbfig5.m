
addpath /mnt/storage/prolatespheroid/
addpath /mnt/storage/modeling_code_current/

% monotype aggregation, aspect ratio=5
gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/2tb_hd/ECHO_STAT_add/20121227_pdfmodel_prosph_ar5_w_respenv_5e4';
fname = 'tx_sys';
tx_opt = 3;  % actual tx
taper_opt = 0;  % no additional taper
N = [1000:100:1400,1600:100:1900,2200,2500,3000];
r = 1;  % scattering strength ratio, 1 for monotype aggregation
ar = 5;  % aspect ratio for prolate spheroid, [] for Rayleigh scatterer
for iN=1:length(N)
    model_n_scat_freqdep_bp_mixed(N(iN),r,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,1/ar);
end


%{

% monotype aggregation, aspect ratio=2
gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121227_pdfmodel_prosph_ar2_w_respenv_5e4';
fname = 'tx_sys';
tx_opt = 3;  % actual tx
taper_opt = 0;  % no additional taper
%N = [10,50,300,500:100:1000];
%N = [100,200,1100:100:1500,2000];
%N = [1600:100:1900];
N = [2200,2500,3000];
r = 1;  % scattering strength ratio, 1 for monotype aggregation
ar = 2;  % aspect ratio for prolate spheroid, [] for Rayleigh scatterer
for iN=1:length(N)
    model_n_scat_freqdep_bp_mixed(N(iN),r,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,1/ar);
end

% monotype aggregation, aspect ratio=10
gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121227_pdfmodel_prosph_ar10_w_respenv_5e4';
fname = 'tx_sys';
tx_opt = 3;  % actual tx
taper_opt = 0;  % no additional taper
%N = [100:100:1500,2000];
%N = [1600:100:1900];
N = [2200,2500,3000];
r = 1;  % scattering strength ratio, 1 for monotype aggregation
ar = 10;  % aspect ratio for prolate spheroid, [] for Rayleigh scatterer
for iN=1:length(N)
    model_n_scat_freqdep_bp_mixed(N(iN),r,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,1/ar);
end

% monotype aggregation, aspect ratio=3
gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121227_pdfmodel_prosph_ar3_w_respenv_5e4';
fname = 'tx_sys';
tx_opt = 3;  % actual tx
taper_opt = 0;  % no additional taper
N = [100:100:1500,2000];
N = [2200,2500,3000];
%N = [1600:100:1900];
%N = [10,50,300,500:100:1000,20:10:40,60:10:90,100,200,400,1200:200:2000];
r = 1;  % scattering strength ratio, 1 for monotype aggregation
ar = 3;  % aspect ratio for prolate spheroid, [] for Rayleigh scatterer
for iN=1:length(N)
    model_n_scat_freqdep_bp_mixed(N(iN),r,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,1/ar);
end

%}