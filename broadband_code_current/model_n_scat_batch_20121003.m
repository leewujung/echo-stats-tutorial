% call model_n_scat_* in a batch

%{
save_opt = 1;
gateLdist = 0.5;
sampleN = 1e4;
use_ideal = 1;  % ideal transmit
taper_opt = 777;  % very narrowband near 50kHz
save_dir = '/mnt/storage/ECHO_STAT/20120814_bb_1e4smpl_2';
fname = 'sqchirp_50kHz_verynarrow16thlen';
N=[1:4,6:9,11:15];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,...
                            taper_opt,save_opt,save_dir,fname);
end
%}

gateLdist = 0.5;
sampleN = 5e4;
use_ideal = 3;  % actual transmit, with system response
taper_opt = 0;  % no

save_opt = 1;
%save_dir = '/mnt/storage/ECHO_STAT/20121007_bb_5e4smpl';
save_dir = '/mnt/storage/ECHO_STAT/20121018_bb_smpl5e4_prosph';

addpath /mnt/storage/prolatespheroid/
addpath /mnt/storage/modeling_code_current/

fname = 'tx_sys';
N=[1500,2000,30,40,60,70,90,150:100:950];
%N=[25,150,710:10:740,760:10:790,810:10:840,860:10:890,910:10:940,960:10:990];
%N = [30:10:100,200:100:1000,2000];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,...
                            taper_opt,save_opt,save_dir,fname);
end
