% call model_n_scat_* in a batch

sampleN = 1e4;
gateLdist = 0.5;
use_ideal = 1;
taper_opt = 0;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20120529_bp_results';

N=[500:100:1000,1500,2000,5000];
fixfreq = 30e3;
for iN=1:length(N)
fixfreq = 30e3;
    model_n_scat_fixed_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt, ...
                          save_opt,save_dir,fixfreq);
fixfreq = 50e3;
    model_n_scat_fixed_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt, ...
                          save_opt,save_dir,fixfreq);
fixfreq = 70e3;
    model_n_scat_fixed_bp(N(iN),sampleN,gateLdist,use_ideal,taper_opt, ...
                          save_opt,save_dir,fixfreq);
end


% 20120223_bp_results: ideal UNtapered tx, freq_dependent bp
% 20120223_bp_results_taper: ideal tapered tx, freq_dependent bp
% 20120511_bp_results: actual tx + system resp, freq-dependent bp
% 20120521_bp_results: actual tx + system resp, fixed bp, 50kHz
% 20120522_bp_results: ideal UNtapered tx, fixed bp, 50kHz
