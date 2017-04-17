% Compare echo pdf with and without added noise

SAVE_DIR = '/mnt/storage/bb_echopdf_figs/echopdfmodel_with_without_noise_20121211';
wnoise_N100to1000_f = '/echopdfmodel_wnoise_N100to1000.mat';

BB_MODEL_DIR = '/mnt/storage/ECHO_STAT/20121116_pdfmodel_w_respenv_5e4';
fpost = 'gateLen0.5_freqDepBP.mat';
N =[10:10:90,100:50:1000];
npt = 1e2;
[pm,pmx,N,s_all] = get_model_curves_kde(BB_MODEL_DIR,N,npt,fpost);
save([SAVE_DIR,'/echopdfmodel_nonoise.mat'],'pm','pmx','N','s_all');

% generate noise-added echo pdf
ds.cruise = 'EN453';
%ds.DATA_DIR = '/mnt/storage/EN453_fish_metadata';
ds.DATA_DIR = '/mnt/storage/EN453_fish_metadata';
ds.DATA_FILE = '20080908_WB1_tri2_020-022';
%ds.DATA_FILE = '20080908_WB1_tri2_023-025';
ds.BP_DIR = '/mnt/storage/calibration_2008/BeamPattern';
ds.CAL_DIR = '/mnt/storage/calibration_2008/cal_result_20120523';
ds.SAVE_DIR = SAVE_DIR;

%addpath(genpath(ds.CODE_DIR))
addpath /mnt/storage/analysis_code_current
addpath /mnt/storage/broadband_code_current
addpath /mnt/storage/net_tows
addpath /mnt/storage/fish_scat_model

% fish box
[d,p] = setBoxDP;   % load in box depth and ping range
% noise box, depth not adjusted
np = [3280,3315];
nd = [30,35];

% Get noise pdf
% not normalized noise pdf
[nrms,nx,ndens,ntar,h_noise] = getNoiseAmp(ds,'AirMarLow',np(1),np(2),nd(1),nd(2),1e2,1);
print([SAVE_DIR,'/noisepdf.png'],'-dpng');
saveas(h_noise,[SAVE_DIR,'/noisepdf.fig'],'fig');

% Echo pdf
[pf,pfx,pfbw,frms,h] = getFishEchoPdf(ds,p(1,1),p(1,2),d(1,1),d(1,2),1e2);
% make noise-added echo pdf model
BB_MODEL_DIR = '/mnt/storage/ECHO_STAT/20121116_pdfmodel_w_respenv_5e4';
fpost = 'gateLen0.5_freqDepBP.mat';
N =[10:10:90,100:50:1000];
npt = 1e2;
[pm,pmx,N,s_all] = get_model_curves_kde(BB_MODEL_DIR,N,npt,fpost,[],frms/nrms);
save([ds.SAVE_DIR,'/echopdfmodel_wnoise_N10-1000.mat'],'pm','pmx','N','s_all','ds');
wnoise_N10to1000_f = 'echopdfmodel_wnoise_N10-1000.mat';

nonoise_f = 'echopdfmodel_nonoise.mat';


% Plot frms nrms comparison
% noise pdf normalized by frms
[ndens_norm,nx_norm,~] = findEchoDist_kde(ntar/frms,1e2);
nvar = var(ntar/frms);  % normalized by frms
nsigma = sqrt(2*nvar/(4-pi));
nrayl_norm = raylpdf(nx_norm,nsigma);
raylx = logspace(-4,2,100);
rayl = raylpdf(raylx,1/sqrt(2));
nrayl_rms = raylpdf(raylx,1/sqrt(2)/(frms/nrms));

h_frms_cmp = figure;
loglog(raylx,rayl,'color',[1 1 1]*180/255,'linewidth',1.5);
hold on
loglog(pfx(1:3:end),pf(1:3:end),'rx','markersize',8);
loglog(nx_norm,ndens_norm,'bx','markersize',8);
loglog(nx_norm,nrayl_norm,'b-','linewidth',1.5);
loglog(raylx,nrayl_rms,'k-','linewidth',1.5);
axis([1e-4 20 1e-5 1e2]);
xlabel('Normazlied echo amplitude','fontsize',12);
ylabel('PDF','fontsize',12);
set(gca,'fontsize',10);
print([SAVE_DIR,'/noisepdf_frms_cmp.png'],'-dpng');
saveas(h_frms_cmp,[SAVE_DIR,'/noisepdf_frms_cmp.fig'],'fig');


% Plot echo pdf comparison
wnoise_N10to1000 = load([SAVE_DIR,'/',wnoise_N10to1000_f]);
nonoise = load([SAVE_DIR,'/',nonoise_f]);

raylx = logspace(-4,2,100);
rayl = raylpdf(raylx,1/sqrt(2));

% N=10;
figure;
loglog(raylx,rayl,'color',[1 1 1]*180/250,'linewidth',1.5);
hold on
loglog(nonoise.pmx,nonoise.pm(:,1),'k','linewidth',0.5);
loglog(wnoise_N10to1000.pmx,wnoise_N10to1000.pm(:,1),'k--','linewidth',3);
xlabel('Normazlied echo amplitude','fontsize',12);
ylabel('PDF','fontsize',12);
set(gca,'fontsize',10);
axis([5e-3 50 1e-5 1e2]);
legend('Rayleigh','No noise','With noise');
title('N=10','fontsize',12)
saveas(gcf,[SAVE_DIR,'/N10_noise_cmp.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/N10_noise_cmp.png'],'png');

% N=50;
figure;
loglog(raylx,rayl,'color',[1 1 1]*180/250,'linewidth',1.5);
hold on
loglog(nonoise.pmx,nonoise.pm(:,5),'k','linewidth',0.5);
loglog(wnoise_N10to1000.pmx,wnoise_N10to1000.pm(:,5),'k--','linewidth',3);
xlabel('Normazlied echo amplitude','fontsize',12);
ylabel('PDF','fontsize',12);
set(gca,'fontsize',10);
axis([5e-3 50 1e-5 1e2]);
legend('Rayleigh','No noise','With noise');
title('N=50','fontsize',12)
saveas(gcf,[SAVE_DIR,'/N50_noise_cmp.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/N50_noise_cmp.png'],'png');

% N=100;
figure;
loglog(raylx,rayl,'color',[1 1 1]*180/250,'linewidth',1.5);
hold on
loglog(nonoise.pmx,nonoise.pm(:,10),'k','linewidth',0.5);
loglog(wnoise_N10to1000.pmx,wnoise_N10to1000.pm(:,10),'k--','linewidth',3);
xlabel('Normazlied echo amplitude','fontsize',12);
ylabel('PDF','fontsize',12);
set(gca,'fontsize',10);
axis([5e-3 50 1e-5 1e2]);
legend('Rayleigh','No noise','With noise');
title('N=100','fontsize',12)
saveas(gcf,[SAVE_DIR,'/N100_noise_cmp.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/N100_noise_cmp.png'],'png');

% N=300;
figure;
loglog(raylx,rayl,'color',[1 1 1]*180/250,'linewidth',1.5);
hold on
loglog(nonoise.pmx,nonoise.pm(:,14),'k','linewidth',0.5);
loglog(wnoise_N10to1000.pmx,wnoise_N10to1000.pm(:,10),'k--','linewidth',3);
xlabel('Normazlied echo amplitude','fontsize',12);
ylabel('PDF','fontsize',12);
set(gca,'fontsize',10);
axis([5e-3 50 1e-5 1e2]);
legend('Rayleigh','No noise','With noise');
title('N=300','fontsize',12)
saveas(gcf,[SAVE_DIR,'/N300_noise_cmp.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/N300_noise_cmp.png'],'png');



