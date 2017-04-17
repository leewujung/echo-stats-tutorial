% Fig: with and without beampattern effects comparison

addpath /mnt/storage/analysis_code_current

DATA_DIR = '/mnt/storage/ECHO_STAT/20121129_with_without_bp_bbmodel';
SAVE_DIR = '/mnt/storage/bb_echopdf_figs/bp_nobp_cmp';

rayl = raylpdf(logspace(-4,2,100),1/sqrt(2));

% N=10
nobp = load([DATA_DIR,'/','tx_nosys_N10_sampleN50000_gateLen0.5_noBP.mat'],'s');
bp = load([DATA_DIR,'/','tx_nosys_N10_sampleN100000_gateLdist0.5_freqDepBP.mat'],'s');
bp.s = bp.s(1:5e4);
[p_nobp10,x_nobp10,~] = findEchoDist_kde(nobp.s/sqrt(mean(nobp.s.^2)),1e2);
[p_bp10,x_bp10,~] = findEchoDist_kde(bp.s/sqrt(mean(bp.s.^2)),1e2);

figure;
loglog(logspace(-4,2,100),rayl,'color',[1 1 1]*180/255,'linewidth',1);
hold on
loglog(x_nobp10,p_nobp10,'k','linewidth',1);
loglog(x_bp10,p_bp10,'k--','linewidth',1);
axis([1e-3 50 1e-5 1e2]);
legend('Rayleigh','No bp','With bp');
title('N=10','fontsize',12);
set(gca,'fontsize',12,'ticklength',[0.02 0.02]);
xlabel('Normalized echo amplitude','fontsize',12);
ylabel('PDF','fontsize',12);
saveas(gcf,[SAVE_DIR,'/bp_nobp_cmp_N10.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/bp_nobp_cmp_N10.png'],'png');


% N=50
nobp = load([DATA_DIR,'/','tx_nosys_N50_sampleN50000_gateLen0.5_noBP.mat'],'s');
bp = load([DATA_DIR,'/','tx_nosys_N50_sampleN500000_gateLdist0.5_freqDepBP.mat'],'s');
bp.s = bp.s(1:5e4);
[p_nobp10,x_nobp10,~] = findEchoDist_kde(nobp.s/sqrt(mean(nobp.s.^2)),1e2);
[p_bp10,x_bp10,~] = findEchoDist_kde(bp.s/sqrt(mean(bp.s.^2)),1e2);

figure;
loglog(logspace(-4,2,100),rayl,'color',[1 1 1]*180/255,'linewidth',1);
hold on
loglog(x_nobp10,p_nobp10,'k','linewidth',1);
loglog(x_bp10,p_bp10,'k--','linewidth',1);
axis([1e-3 50 1e-5 1e2]);
legend('Rayleigh','No bp','With bp');
title('N=50','fontsize',12);
set(gca,'fontsize',12,'ticklength',[0.02 0.02]);
xlabel('Normalized echo amplitude','fontsize',12);
ylabel('PDF','fontsize',12);
saveas(gcf,[SAVE_DIR,'/bp_nobp_cmp_N50.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/bp_nobp_cmp_N50.png'],'png');


% N=300
nobp = load([DATA_DIR,'/','tx_nosys_N300_sampleN50000_gateLen0.5_noBP.mat'],'s');
bp = load([DATA_DIR,'/','tx_nosys_N300_sampleN100000_gateLdist0.5_freqDepBP.mat'],'s');
bp.s = bp.s(1:5e4);
[p_nobp10,x_nobp10,~] = findEchoDist_kde(nobp.s/sqrt(mean(nobp.s.^2)),1e2);
[p_bp10,x_bp10,~] = findEchoDist_kde(bp.s/sqrt(mean(bp.s.^2)),1e2);

figure;
loglog(logspace(-4,2,100),rayl,'color',[1 1 1]*180/255,'linewidth',1);
hold on
loglog(x_nobp10,p_nobp10,'k','linewidth',1);
loglog(x_bp10,p_bp10,'k--','linewidth',1);
axis([1e-3 50 1e-5 1e2]);
legend('Rayleigh','No bp','With bp');
title('N=300','fontsize',12);
set(gca,'fontsize',12,'ticklength',[0.02 0.02]);
xlabel('Normalized echo amplitude','fontsize',12);
ylabel('PDF','fontsize',12);
saveas(gcf,[SAVE_DIR,'/bp_nobp_cmp_N300.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/bp_nobp_cmp_N300.png'],'png');


