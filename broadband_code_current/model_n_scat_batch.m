% call model_n_scat_* in a batch

sampleN = 1e4;
gateLdist = 0.5;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20120731_bp_results';

%{
taper_opt = 2;  % taper, only HF part exists
use_ideal = 1;  % square chirp
fname = 'hann_hf_taper';
N=[10,50,100,200,300,500];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt,save_opt,save_dir,fname);
end
%}

%{
taper_opt = 3;  % taper, only LF part exists
use_ideal = 1;  % square chirp
fname = 'hann_lf_taper';
N=[10,50,100,200,300,500];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt,save_opt,save_dir,fname);
end
%}

%{
taper_opt = 0;  % no taper
use_ideal = 1;  % square chirp
fname = 'no_taper';
N=[10,50,100,200,300,500];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt,save_opt,save_dir,fname);
end
%}


%{
taper_opt = 0;  % no taper
use_ideal = 3;  % use actual tx + sys resp
fname = 'actual_w_sys_no_taper';
N=[10,50,100,200,300,500];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt,save_opt,save_dir,fname);
end

taper_opt = 0;  % no taper
use_ideal = 2;  % use ideal tx, no sys resp
fname = 'ideal_no_sys_no_taper';
N=[10,50,100,200,300,500];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt,save_opt,save_dir,fname);
end

taper_opt = 1;  % full Gaussian taper
use_ideal = 1;  % square chirp
fname = 'square_chirp_gaussian_taper';
N=[10,50,100,200,300,500];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt,save_opt,save_dir,fname);
end
%}

%{
disp('middle Hann taper, freq-dep bp');
taper_opt = 4;  % middle Hann taper
use_ideal = 1;  % square chirp
fname = 'middle_hann_taper';
N=[50,200,500];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt,save_opt,save_dir,fname);
end

disp('upper Hann taper, freq-dep bp');
taper_opt = 5;  % upper Hann taper
use_ideal = 1;  % square chirp
fname = 'upper_hann_taper';
N=[50,200,500];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt,save_opt,save_dir,fname);
end

disp('lower Hann taper, freq-dep bp');
taper_opt = 6;  % lower Hann taper
use_ideal = 1;  % square chirp
fname = 'lower_hann_taper';
N=[50,200,500];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt,save_opt,save_dir,fname);
end
%}


%{
disp('sq, fixed freq');
taper_opt = 0;  % no taper
use_ideal = 1;  % square chirp
fname = 'square_chirp_no_taper';
N=[50,200,500];
for iN=1:length(N)
    model_n_scat_fixed_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt, ...
                            save_opt,save_dir,30e3,fname);
    model_n_scat_fixed_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt, ...
                            save_opt,save_dir,70e3,fname);
    model_n_scat_fixed_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt, ...
                            save_opt,save_dir,30e3,fname);
end
%}

%{
disp('mid narrow Hann taper, freq-dep bp');
taper_opt = 7;  % lower Hann taper
use_ideal = 1;  % square chirp
fname = 'middle_hann_taper_narrow';
N=[10,100,300];
%N=[50,200,500];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt,save_opt,save_dir,fname);
end

disp('mid wide Hann taper, freq-dep bp');
taper_opt = 9;  % lower Hann taper
use_ideal = 1;  % square chirp
fname = 'middle_hann_taper_wide';
N=[10,100,300];
%N=[50,200,500];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt,save_opt,save_dir,fname);
end

disp('mid narrow Hann taper, freq-dep bp');
taper_opt = 7;  % lower Hann taper
use_ideal = 1;  % square chirp
fname = 'middle_hann_taper_narrow';
N=[50,200,500];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt,save_opt,save_dir,fname);
end

disp('mid wide Hann taper, freq-dep bp');
taper_opt = 9;  % lower Hann taper
use_ideal = 1;  % square chirp
fname = 'middle_hann_taper_wide';
N=[50,200,500];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt,save_opt,save_dir,fname);
end
%}



disp('sq chirp, no bp');
taper_opt = 0;  % no taper
use_ideal = 1;  % square chirp
fname = 'square_chirp_no_bp';
N=[10,100,300,50,200,500];
for iN=1:length(N)
    model_n_scat_no_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt, ...
                            save_opt,save_dir,fname);
end
