% Test the influence of the angle of orientation distribution of fish

N = [10,20,50,100,200,500];
mix_r = 1;
ns = 5e4;
glen = 0.5;
save_opt = 1;
taper = 0;  % taper transmit signal
tx = 3;
sdir = '/mnt/storage/ECHO_STAT/20130802_test_pdfs';
figdir = '/mnt/storage/ECHO_STAT/20130805_fish_distr_effect';

for iN=1:length(N)

%{
% Rayleigh
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'rayleigh',0);

% Fish, 22cm only
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std5',...
         1,[0.22,1],[0,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean-5_std5',...
         1,[0.22,1],[-5,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean-10_std5',...
         1,[0.22,1],[-10,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std10',...
         1,[0.22,1],[0,10]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std15',...
         1,[0.22,1],[0,15]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean0_std20',...
         1,[0.22,1],[0,20]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean5_std0',...
         1,[0.22,1],[5,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean10_std0',...
         1,[0.22,1],[10,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_len22_mean20_std0',...
         1,[0.22,1],[20,0]);
%}

% Fish, empirincal distribution
%{
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean0_std5',...
         1,[],[0,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean-5_std5',...
         1,[],[-5,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean-10_std5',...
         1,[],[-10,5]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean0_std10',...
         1,[],[0,10]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean0_std15',...
         1,[],[0,15]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean0_std20',...
         1,[],[0,20]);

model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean5_std0',...
         1,[],[5,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean10_std0',...
         1,[],[10,0]);
model_n_scat_freqdep_bp_mixed_freqdomain_wfish(N(iN),mix_r,ns, ...
         glen,tx,taper,save_opt,sdir,'fish_lendistr_mean20_std0',...
         1,[],[20,0]);
%}



stail = ['_sampleN',num2str(ns),'_gateLen0.5_freqDepBP.mat'];
cd('/mnt/storage/ECHO_STAT/20130802_test_pdfs');

R = load(['rayleigh_N',num2str(N(iN)),stail]);

F22_0_5  = load(['fish_len22_mean0_std5_N',num2str(N(iN)),stail]);
F22_0_10 = load(['fish_len22_mean0_std10_N',num2str(N(iN)),stail]);
F22_0_15 = load(['fish_len22_mean0_std15_N',num2str(N(iN)),stail]);
F22_0_20 = load(['fish_len22_mean0_std20_N',num2str(N(iN)),stail]);
F22_5_5  = load(['fish_len22_mean-5_std5_N',num2str(N(iN)),stail]);
F22_10_5 = load(['fish_len22_mean-10_std5_N',num2str(N(iN)),stail]);

F_0_5  = load(['fish_lendistr_mean0_std5_N',num2str(N(iN)),stail]);
F_0_10 = load(['fish_lendistr_mean0_std10_N',num2str(N(iN)),stail]);
F_0_15 = load(['fish_lendistr_mean0_std15_N',num2str(N(iN)),stail]);
F_0_20 = load(['fish_lendistr_mean0_std20_N',num2str(N(iN)),stail]);
F_5_5  = load(['fish_lendistr_mean-5_std5_N',num2str(N(iN)),stail]);
F_10_5 = load(['fish_lendistr_mean-10_std5_N',num2str(N(iN)),stail]);

F22_5_0  = load(['fish_len22_mean5_std0_N',num2str(N(iN)),stail]);
F22_10_0 = load(['fish_len22_mean10_std0_N',num2str(N(iN)),stail]);
F22_20_0 = load(['fish_len22_mean20_std0_N',num2str(N(iN)),stail]);

F_5_0  = load(['fish_lendistr_mean5_std0_N',num2str(N(iN)),stail]);
F_10_0 = load(['fish_lendistr_mean10_std0_N',num2str(N(iN)),stail]);
F_20_0 = load(['fish_lendistr_mean20_std0_N',num2str(N(iN)),stail]);

raylx = logspace(-5,50,1000);
rayl = raylpdf(raylx,1/sqrt(2));

[pmR,pmxR] = findEchoDist_kde(R.s/sqrt(mean(R.s.^2)),200);

[pmF22_0_5,pmxF22_0_5] = findEchoDist_kde(F22_0_5.s/sqrt(mean(F22_0_5.s.^2)),200);
[pmF22_0_10,pmxF22_0_10] = findEchoDist_kde(F22_0_10.s/sqrt(mean(F22_0_10.s.^2)),200);
[pmF22_0_15,pmxF22_0_15] = findEchoDist_kde(F22_0_15.s/sqrt(mean(F22_0_15.s.^2)),200);
[pmF22_0_20,pmxF22_0_20] = findEchoDist_kde(F22_0_20.s/sqrt(mean(F22_0_20.s.^2)),200);
[pmF22_5_5,pmxF22_5_5] = findEchoDist_kde(F22_5_5.s/sqrt(mean(F22_5_5.s.^2)),200);
[pmF22_10_5,pmxF22_10_5] = findEchoDist_kde(F22_10_5.s/sqrt(mean(F22_10_5.s.^2)),200);

[pmF22_5_0,pmxF22_5_0] = findEchoDist_kde(F22_5_0.s/sqrt(mean(F22_5_0.s.^2)),200);
[pmF22_10_0,pmxF22_10_0] = findEchoDist_kde(F22_10_0.s/sqrt(mean(F22_10_0.s.^2)),200);
[pmF22_20_0,pmxF22_20_0] = findEchoDist_kde(F22_20_0.s/sqrt(mean(F22_20_0.s.^2)),200);

[pmF_0_5,pmxF_0_5] = findEchoDist_kde(F_0_5.s/sqrt(mean(F_0_5.s.^2)),200);
[pmF_0_10,pmxF_0_10] = findEchoDist_kde(F_0_10.s/sqrt(mean(F_0_10.s.^2)),200);
[pmF_0_15,pmxF_0_15] = findEchoDist_kde(F_0_15.s/sqrt(mean(F_0_15.s.^2)),200);
[pmF_0_20,pmxF_0_20] = findEchoDist_kde(F_0_20.s/sqrt(mean(F_0_20.s.^2)),200);
[pmF_5_5,pmxF_5_5] = findEchoDist_kde(F_5_5.s/sqrt(mean(F_5_5.s.^2)),200);
[pmF_10_5,pmxF_10_5] = findEchoDist_kde(F_10_5.s/sqrt(mean(F_10_5.s.^2)),200);

[pmF_5_0,pmxF_5_0] = findEchoDist_kde(F_5_0.s/sqrt(mean(F_5_0.s.^2)),200);
[pmF_10_0,pmxF_10_0] = findEchoDist_kde(F_10_0.s/sqrt(mean(F_10_0.s.^2)),200);
[pmF_20_0,pmxF_20_0] = findEchoDist_kde(F_20_0.s/sqrt(mean(F_20_0.s.^2)),200);


%% Plot
% No fish length variation, only angle of orientation distribution
figure;
loglog(raylx,rayl,'k--');
hold on
loglog(pmxR,pmR,'k');
loglog(pmxF22_0_5,pmF22_0_5,'b');
loglog(pmxF22_0_10,pmF22_0_10,'r');
loglog(pmxF22_0_15,pmF22_0_15,'color',[0,77,0]/255);
loglog(pmxF22_0_20,pmF22_0_20,'color',[204,153,0]/255);
legend('Rayleigh distr','Rayleigh scatterer',...
       '22cm fish, mean 0, std 5',...
       '22cm fish, mean 0, std 10',...
       '22cm fish, mean 0, std 15',...
       '22cm fish, mean 0, std 20');
axis([1e-3 1e2 1e-7 1e2])
title('Fish 22cm, same mean, diff std');
saveas(gcf,[figdir,'/fish22cm_sameMean_diffStd_N',num2str(N(iN)),'.png'],'png');

figure;
loglog(raylx,rayl,'k');
hold on
loglog(pmxR,pmR,'k');
loglog(pmxF22_0_5,pmF22_0_5,'b');
loglog(pmxF22_5_5,pmF22_5_5,'r');
loglog(pmxF22_10_5,pmF22_10_5,'color',[0,77,0]/255);
legend('Rayleigh distr','Rayleigh scatterer',...
       '22cm fish, mean 0, std 5',...
       '22cm fish, mean 5, std 5',...
       '22cm fish, mean 10, std 5');
axis([1e-3 1e2 1e-7 1e2])
title('Fish 22cm, same std, diff mean');
saveas(gcf,[figdir,'/fish22cm_sameStd_diffMean_N',num2str(N(iN)),'.png'],'png');

% Both angle of orientation and fish length distribution
figure;
loglog(raylx,rayl,'k--');
hold on
loglog(pmxR,pmR,'k');
loglog(pmxF_0_5,pmF_0_5,'b');
loglog(pmxF_0_10,pmF_0_10,'r');
loglog(pmxF_0_15,pmF_0_15,'color',[0,77,0]/255);
loglog(pmxF_0_20,pmF_0_20,'color',[204,153,0]/255);
legend('Rayleigh distr','Rayleigh scatterer',...
       'fish all length, mean 0, std 5',...
       'fish all length, mean 0, std 10',...
       'fish all length, mean 0, std 15',...
       'fish all length, mean 0, std 20');
axis([1e-3 1e2 1e-7 1e2])
title('Fish all length, same mean, diff std');
saveas(gcf,[figdir,'/fishAllLen_sameMean_diffStd_N',num2str(N(iN)),'.png'],'png');

figure;
loglog(raylx,rayl,'k');
hold on
loglog(pmxR,pmR,'k');
loglog(pmxF_0_5,pmF_0_5,'b');
loglog(pmxF_5_5,pmF_5_5,'r');
loglog(pmxF_10_5,pmF_10_5,'color',[0,77,0]/255);
legend('Rayleigh distr','Rayleigh scatterer',...
       'fish all length, mean 0, std 5',...
       'fish all length, mean 5, std 5',...
       'fish all length, mean 10, std 5');
axis([1e-3 1e2 1e-7 1e2])
title('Fish all length, same std, diff mean');
saveas(gcf,[figdir,'/fishAllLen_sameStd_diffMean_N',num2str(N(iN)),'.png'],'png');

% No angle of orientation variation, no fish length distribution
figure;
loglog(raylx,rayl,'k');
hold on
loglog(pmxR,pmR,'k');
loglog(pmxF_5_0,pmF_5_0,'b');
loglog(pmxF_10_0,pmF_10_0,'r');
loglog(pmxF_20_0,pmF_20_0,'color',[0,77,0]/255);
legend('Rayleigh distr','Rayleigh scatterer',...
       'fish all length, mean 0, std 0',...
       'fish all length, mean 10, std 0',...
       'fish all length, mean 20, std 0');
axis([1e-3 1e2 1e-7 1e2])
title('Fish all length, same angle only');
saveas(gcf,[figdir,'/fishAllLen_sameAngle_N',num2str(N(iN)),'.png'],'png');

% No angle of orientation variation, only fish length distribution
figure;
loglog(raylx,rayl,'k');
hold on
loglog(pmxR,pmR,'k');
loglog(pmxF22_5_0,pmF22_5_0,'b');
loglog(pmxF22_10_0,pmF22_10_0,'r');
loglog(pmxF22_20_0,pmF22_20_0,'color',[0,77,0]/255);
legend('Rayleigh distr','Rayleigh scatterer',...
       '22cm fish, mean 0, std 0',...
       '22cm fish, mean 10, std 0',...
       '22cm fish, mean 20, std 0');
axis([1e-3 1e2 1e-7 1e2])
title('22cm fish, same angle only');
saveas(gcf,[figdir,'/fish22cm_sameAngle_N',num2str(N(iN)),'.png'],'png');

% No angle of orientation distribution, fish length distribution comparison
figure;
loglog(raylx,rayl,'k');
hold on
loglog(pmxR,pmR,'k');
loglog(pmxF22_5_0,pmF22_5_0,'b');
loglog(pmxF_5_0,pmF_5_0,'r');
legend('Rayleigh distr','Rayleigh scatterer',...
       '22cm fish, mean 5, std 0',...
       'fish all length, mean 5, std 0');
axis([1e-3 1e2 1e-7 1e2])
title('same angle, diff fish length distr');
saveas(gcf,[figdir,'/fishLenDistrCmp_mean5_std0_N',num2str(N(iN)),'.png'],'png');

figure;
loglog(raylx,rayl,'k');
hold on
loglog(pmxR,pmR,'k');
loglog(pmxF22_10_0,pmF22_10_0,'b');
loglog(pmxF_10_0,pmF_10_0,'r');
legend('Rayleigh distr','Rayleigh scatterer',...
       '22cm fish, mean 10, std 0',...
       'fish all length, mean 10, std 0');
axis([1e-3 1e2 1e-7 1e2])
title('same angle, diff fish length distr');
saveas(gcf,[figdir,'/fishLenDistrCmp_mean10_std0_N',num2str(N(iN)),'.png'],'png');

figure;
loglog(raylx,rayl,'k');
hold on
loglog(pmxR,pmR,'k');
loglog(pmxF22_20_0,pmF22_20_0,'b');
loglog(pmxF_20_0,pmF_20_0,'r');
legend('Rayleigh distr','Rayleigh scatterer',...
       '22cm fish, mean 20, std 0',...
       'fish all length, mean 20, std 0');
axis([1e-3 1e2 1e-7 1e2])
title('same angle, diff fish length distr');
saveas(gcf,[figdir,'/fishLenDistrCmp_mean20_std0_N',num2str(N(iN)),'.png'],'png');


end