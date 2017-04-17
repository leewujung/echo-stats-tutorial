% 2012 11 07  Plot bpir and assoc. xcorr results using different signals


SAVE_DIR = '/mnt/storage/bb_echopdf_figs/bpir';

%% Load beampattern ir
bp_file = '/mnt/storage/ECHO_STAT/20120223_bp_results/bpir_a0.054m_dtheta0.001pi_fmax10000kHz.mat';
BPIR = load(bp_file);
dt = diff(BPIR.t_all(1:2));
bpirL = length(BPIR.t_all);
bpirHalfL = (length(BPIR.t_all)+1)/2;  % half length of the bpir

%% Load signals
[y,t_y] = gen_tx(1,dt);
[p,q] = rat(diff(t_y(1:2))/dt);
y_sqchirp = resample(y,p,q);
t_y_sqchirp = ((1:length(y))-1)*dt;

[y,t_y] = gen_tx(2,dt);
[p,q] = rat(diff(t_y(1:2))/dt);
y_idealtx = resample(y,p,q);
t_y_idealtx = ((1:length(y_idealtx))-1)*dt;

[y,t_y] = gen_tx(3,dt);
[p,q] = rat(diff(t_y(1:2))/dt);
y_actualtx = resample(y,p,q);
t_y_actualtx = ((1:length(y_actualtx))-1)*dt;

% window the chirp
win = win_chirp(11,y_sqchirp);
y_hihann = y_sqchirp.*win;
win = win_chirp(12,y_sqchirp);
y_lohann = y_sqchirp.*win;
win = win_chirp(7,y_sqchirp);
y_midhnarrow = y_sqchirp.*win;
win = win_chirp(9,y_sqchirp);
y_midhwide = y_sqchirp.*win;


len = length(y_idealtx)-length(y_sqchirp);
y_sqchirp = [y_sqchirp,zeros(1,len)];
y_hihann = [y_hihann,zeros(1,len)];
y_lohann = [y_lohann,zeros(1,len)];
y_midhnarrow = [y_midhnarrow,zeros(1,len)];
y_midhwide = [y_midhwide,zeros(1,len)];

t_y = t_y_idealtx;

[sig_xcorr(1,:),t_ya] = xcorr(y_sqchirp);
sig_xcorr(2,:) = xcorr(y_idealtx);
sig_xcorr(3,:) = xcorr(y_actualtx);
sig_xcorr(4,:) = xcorr(y_hihann);
sig_xcorr(5,:) = xcorr(y_lohann);
sig_xcorr(6,:) = xcorr(y_midhnarrow);
sig_xcorr(7,:) = xcorr(y_midhwide);


% Select bpir

theta_sel = [5,10,30]/180*pi;
for iT = 1:length(theta_sel)
    [~,idx] = min(abs(BPIR.theta-theta_sel(iT)));
    theta_idx(iT) = idx;
    bpir(iT,:) = BPIR.bpir_all(idx,:);
    for iM=1:7
        [xya(iT,iM,:),t_xya] = xcorr(sig_xcorr(iM,:),bpir(iT,:));
    end
end

for iT = 1:length(theta_sel)
    xya_env(iT,:,:) = abs(hilbert(squeeze(xya(iT,:,:)).')).';
end
m = max(squeeze(xya_env(1,:,:)),[],2);

[~,idx] = max(BPIR.bpir_all(1,:));
mididx = ceil((length(t_y)+1)/2)*2+idx;
idx = mididx+(-5000:5000);
xya_env = xya_env(:,:,idx);

xya_env = xya_env./permute(repmat(m,[1,3,length(idx)]),[2,1,3]);

t_bpir = BPIR.t_all*1e3;
t_xya = (-5000:5000)*dt*1e3;



% bpir
figure;
plot(t_bpir,BPIR.bpir_all(theta_idx(1),:),'k','linewidth',1);
hold on
plot(t_bpir,BPIR.bpir_all(theta_idx(2),:),'k--','linewidth',1);
plot(t_bpir,BPIR.bpir_all(theta_idx(3),:),'k','linewidth',2);
legend('\theta=5^{\circ}','\theta=10^{\circ}','\theta=30^{\circ}')
axis([-0.04 0.04 0 5e5])
set(gca,'xtick',-0.04:0.02:0.04,'ytick',(0:5)*1e5);
set(gca,'yticklabel',num2str((0:5)'),'fontsize',12);
xlabel('Time (ms)','fontsize',14);
ylabel('Arbitrary relative scale','fontsize',14);
set(gca,'TickLength',[0.02,0]);
set(gca,'Layer','Top');
saveas(gcf,[SAVE_DIR,'/bp.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/bp.png'],'png');

% square chirp
figure;
plot(t_xya,squeeze(xya_env(1,1,:)),'k','linewidth',1);
hold on
plot(t_xya,squeeze(xya_env(2,1,:)),'k--','linewidth',1);
plot(t_xya,squeeze(xya_env(3,1,:)),'k','linewidth',2);
legend('\theta=5^{\circ}','\theta=10^{\circ}','\theta=30^{\circ}');
axis([-0.15 0.15 0 1.1])
set(gca,'xtick',-0.15:0.05:0.15,'ytick',0:0.5:1);
set(gca,'yticklabel',num2str((0:0.5:1)'),'fontsize',12);
xlabel('Time (ms)','fontsize',14);
ylabel('Normazlied envelope','fontsize',14);
set(gca,'TickLength',[0.02,0]);
set(gca,'Layer','Top');
saveas(gcf,[SAVE_DIR,'/bpxcorr_sqchirp.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/bpxcorr_sqchirp.png'],'png');

% ideal transmit signal
figure;
plot(t_xya,squeeze(xya_env(1,2,:)),'k','linewidth',1);
hold on
plot(t_xya,squeeze(xya_env(2,2,:)),'k--','linewidth',1);
plot(t_xya,squeeze(xya_env(3,2,:)),'k','linewidth',2);
legend('\theta=5^{\circ}','\theta=10^{\circ}','\theta=30^{\circ}')
axis([-0.15 0.15 0 1.1])
set(gca,'xtick',-0.15:0.05:0.15,'ytick',0:0.5:1);
set(gca,'yticklabel',num2str((0:0.5:1)'),'fontsize',12);
xlabel('Time (ms)','fontsize',14);
ylabel('Normazlied envelope','fontsize',14);
set(gca,'TickLength',[0.02,0]);
set(gca,'Layer','Top');
saveas(gcf,[SAVE_DIR,'/bpxcorr_idealtx.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/bpxcorr_idealtx.png'],'png');

% actual transmit signal
figure;
plot(t_xya,squeeze(xya_env(1,3,:)),'k','linewidth',1);
hold on
plot(t_xya,squeeze(xya_env(2,3,:)),'k--','linewidth',1);
plot(t_xya,squeeze(xya_env(3,3,:)),'k','linewidth',2);
legend('\theta=5^{\circ}','\theta=10^{\circ}','\theta=30^{\circ}')
axis([-0.15 0.15 0 1.1])
set(gca,'xtick',-0.15:0.05:0.15,'ytick',0:0.5:1);
set(gca,'yticklabel',num2str((0:0.5:1)'),'fontsize',12);
xlabel('Time (ms)','fontsize',14);
ylabel('Normazlied envelope','fontsize',14);
set(gca,'TickLength',[0.02,0]);
set(gca,'Layer','Top');
saveas(gcf,[SAVE_DIR,'/bpxcorr_actualtx.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/bpxcorr_actualtx.png'],'png');


% high hann win
figure;
plot(t_xya,squeeze(xya_env(1,4,:)),'k','linewidth',1);
hold on
plot(t_xya,squeeze(xya_env(2,4,:)),'k--','linewidth',1);
plot(t_xya,squeeze(xya_env(3,4,:)),'k','linewidth',2);
legend('\theta=5^{\circ}','\theta=10^{\circ}','\theta=30^{\circ}')
axis([-0.15 0.15 0 1.1])
set(gca,'xtick',-0.15:0.05:0.15,'ytick',0:0.5:1);
set(gca,'yticklabel',num2str((0:0.5:1)'),'fontsize',12);
xlabel('Time (ms)','fontsize',14);
ylabel('Normazlied envelope','fontsize',14);
set(gca,'TickLength',[0.02,0]);
set(gca,'Layer','Top');
saveas(gcf,[SAVE_DIR,'/bpxcorr_hihann.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/bpxcorr_hihann.png'],'png');


% low hann win
figure;
plot(t_xya,squeeze(xya_env(1,5,:)),'k','linewidth',1);
hold on
plot(t_xya,squeeze(xya_env(2,5,:)),'k--','linewidth',1);
plot(t_xya,squeeze(xya_env(3,5,:)),'k','linewidth',2);
legend('\theta=5^{\circ}','\theta=10^{\circ}','\theta=30^{\circ}')
axis([-0.15 0.15 0 1.1])
set(gca,'xtick',-0.15:0.05:0.15,'ytick',0:0.5:1);
set(gca,'yticklabel',num2str((0:0.5:1)'),'fontsize',12);
xlabel('Time (ms)','fontsize',14);
ylabel('Normazlied envelope','fontsize',14);
set(gca,'TickLength',[0.02,0]);
set(gca,'Layer','Top');
saveas(gcf,[SAVE_DIR,'/bpxcorr_lohann.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/bpxcorr_lohann.png'],'png');


% mid hann narrow win
figure;
plot(t_xya,squeeze(xya_env(1,6,:)),'k','linewidth',1);
hold on
plot(t_xya,squeeze(xya_env(2,6,:)),'k--','linewidth',1);
plot(t_xya,squeeze(xya_env(3,6,:)),'k','linewidth',2);
legend('\theta=5^{\circ}','\theta=10^{\circ}','\theta=30^{\circ}')
axis([-0.15 0.15 0 1.1])
set(gca,'xtick',-0.15:0.05:0.15,'ytick',0:0.5:1);
set(gca,'yticklabel',num2str((0:0.5:1)'),'fontsize',12);
xlabel('Time (ms)','fontsize',14);
ylabel('Normazlied envelope','fontsize',14);
set(gca,'TickLength',[0.02,0]);
set(gca,'Layer','Top');
saveas(gcf,[SAVE_DIR,'/bpxcorr_midhnarrow.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/bpxcorr_midhnarrow.png'],'png');


% mid hann wide win
figure;
plot(t_xya,squeeze(xya_env(1,7,:)),'k','linewidth',1);
hold on
plot(t_xya,squeeze(xya_env(2,7,:)),'k--','linewidth',1);
plot(t_xya,squeeze(xya_env(3,7,:)),'k','linewidth',2);
legend('\theta=5^{\circ}','\theta=10^{\circ}','\theta=30^{\circ}')
axis([-0.15 0.15 0 1.1])
set(gca,'xtick',-0.15:0.05:0.15,'ytick',0:0.5:1);
set(gca,'yticklabel',num2str((0:0.5:1)'),'fontsize',12);
xlabel('Time (ms)','fontsize',14);
ylabel('Normazlied envelope','fontsize',14);
set(gca,'TickLength',[0.02,0]);
set(gca,'Layer','Top');
saveas(gcf,[SAVE_DIR,'/bpxcorr_midhwide.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/bpxcorr_midhwide.png'],'png');
