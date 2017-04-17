% call model_n_scat_* in a batch


%{
save_dir = '/mnt/storage/ECHO_STAT/20120814_bb_1e4smpl_2';
fname = 'tx_sys';
%N=[10,50,100,200,300,500,20,80,400,600:100:1000,2000];
N=[30,40,50,70,90];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,...
                            taper_opt,save_opt,save_dir,fname);
end
%}

%{
use_ideal = 1;  % ideal signal
taper_opt = 7;  % taper very narrow hann window
save_dir = '/mnt/storage/ECHO_STAT/20120814_bb_1e4smpl_2';
fname = 'sqchirp_narrow1_hann';
N=[50];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,...
                            taper_opt,save_opt,save_dir,fname);
end

use_ideal = 1;  % ideal signal
taper_opt = 10;  % taper very wide hann window
save_dir = '/mnt/storage/ECHO_STAT/20120814_bb_1e4smpl_2';
fname = 'sqchirp_wide2_hann';
N=[50];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,...
                            taper_opt,save_opt,save_dir,fname);
end
%}


%{
sampleN = 4e5;
use_ideal = 2;  % ideal transmit, no system response
taper_opt = 0;  % no
fname = 'tx_nosys';
N=[10];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,...
                            taper_opt,save_opt,save_dir,fname);
end

sampleN = 4e5;
use_ideal = 1;  % square chirp
taper_opt = 11;  % upper Hann wide
fname = 'sqchirp_upperHwide';
N=[10];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,...
                            taper_opt,save_opt,save_dir,fname);
end

sampleN = 4e5;
use_ideal = 1;  % square chirp
taper_opt = 12;  % lower Hann wide
fname = 'sqchirp_lowerHwide';
N=[10];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,...
                            taper_opt,save_opt,save_dir,fname);
end
%}


sampleN = 1e4;
use_ideal = 3;  % actual transmit, with system response
taper_opt = 0;  % no
save_dir = '/mnt/storage/ECHO_STAT/20120814_bb_1e4smpl_2';
fname = 'tx_sys';
N=[5:10:105,150:100:950,1500];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,...
                            taper_opt,save_opt,save_dir,fname);
end
