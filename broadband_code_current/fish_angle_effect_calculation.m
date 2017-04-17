% Test the influence of the angle of orientation distribution of fish

N = [10,20,50,100,200,500];
%N = 10;
mix_r = 1;
ns = 5e4;
glen = 0.5;
save_opt = 1;
taper = 0;  % taper transmit signal
tx = 3;
sdir = '/mnt/storage/ECHO_STAT/20130806_fish_distr_effect_matfiles';

for iN=1:length(N)

% Rayleigh
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'rayleigh',[],0);

% Fish, 22cm only
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std5',...
         [],1,[0.22,1],[0,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean5_std5',...
         [],1,[0.22,1],[5,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean10_std5',...
         [],1,[0.22,1],[10,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean20_std5',...
         [],1,[0.22,1],[20,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean50_std5',...
         [],1,[0.22,1],[50,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std10',...
         [],1,[0.22,1],[0,10]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std15',...
         [],1,[0.22,1],[0,15]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std20',...
         [],1,[0.22,1],[0,20]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean5_std0',...
         [],1,[0.22,1],[5,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean10_std0',...
         [],1,[0.22,1],[10,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean20_std0',...
         [],1,[0.22,1],[20,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean50_std0',...
         [],1,[0.22,1],[50,0]);

%{
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std0_bpa_all',...
         [],1,[0.22,1],[0,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std0_bpa_10',...
         10,1,[0.22,1],[0,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std0_bpa_20',...
         20,1,[0.22,1],[0,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std0_bpa_30',...
         30,1,[0.22,1],[0,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std0_bpa_40',...
         40,1,[0.22,1],[0,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std0_bpa_50',...
         50,1,[0.22,1],[0,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std0_bpa_60',...
         60,1,[0.22,1],[0,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std0_bpa_70',...
         70,1,[0.22,1],[0,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std0_bpa_80',...
         80,1,[0.22,1],[0,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std0_bpa_90',...
         90,1,[0.22,1],[0,0]);
%}

% Fish, empirincal distribution
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean0_std5',...
         [],1,[],[0,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean5_std5',...
         [],1,[],[5,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean10_std5',...
         [],1,[],[10,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean20_std5',...
         [],1,[],[20,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean50_std5',...
         [],1,[],[50,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean0_std10',...
         [],1,[],[0,10]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean0_std15',...
         [],1,[],[0,15]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean0_std20',...
         [],1,[],[0,20]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean5_std0',...
         [],1,[],[5,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean10_std0',...
         [],1,[],[10,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean20_std0',...
         [],1,[],[20,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean50_std0',...
         [],1,[],[50,0]);

end