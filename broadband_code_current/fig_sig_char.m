% 2012 11 08  Plot signal spectral and temporal features


SAVE_DIR = '/mnt/storage/bb_echopdf_figs/sig_char';

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

sig = [y_sqchirp;y_idealtx;y_actualtx;...
       y_hihann;y_lohann;y_midhnarrow;y_midhwide];
sig_fft = fft(sig.').';
sig_fft_max = max(10*log10(abs(sig_fft).^2),[],2);

freq_sig = linspace(0,1/dt/2,ceil(size(sig_fft,2)/2));
fidx = 1:length(freq_sig);
t_y = t_y_idealtx;

%% Xcorr
[sig_xcorr(1,:),t_ya] = xcorr(y_sqchirp);
sig_xcorr(2,:) = xcorr(y_idealtx);
sig_xcorr(3,:) = xcorr(y_actualtx);
sig_xcorr(4,:) = xcorr(y_hihann);
sig_xcorr(5,:) = xcorr(y_lohann);
sig_xcorr(6,:) = xcorr(y_midhnarrow);
sig_xcorr(7,:) = xcorr(y_midhwide);

sig_xcorr_env = abs(hilbert(sig_xcorr.')).';
m = max(sig_xcorr_env,[],2);
sig_xcorr_env = sig_xcorr_env./repmat(m,1,size(sig_xcorr_env,2));
mididx = floor(size(sig_xcorr,2)/2);
idx = mididx+(-5000:5000);
t_xcorr = (-5000:5000)*dt;
sig_xcorr = sig_xcorr(:,idx);
sig_xcorr_env = sig_xcorr_env(:,idx);



%% Plot
% Spectrum
figure;
plot(freq_sig/1e3,10*log10(abs(sig_fft(5,fidx)).^2)-sig_fft_max(5),...
     'k-','linewidth',1);
hold on
plot(freq_sig/1e3,10*log10(abs(sig_fft(4,fidx)).^2)-sig_fft_max(4),...
     'k--','linewidth',3);
axis([20 80 -40 5]);
legend('Low Hann win','High Hann win','location','south');
set(gca,'xtick',20:20:80,'ytick',-40:20:0);
set(gca,'fontsize',12,'ticklength',[1 1]*0.02);
xlabel('Frequency (kHz)','fontsize',14);
ylabel('Arbitrarily normalized unit (dB)','fontsize',14);
saveas(gcf,[SAVE_DIR,'/spectrum_hihann_lohann.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/spectrum_hihann_lohann.png'],'png');

figure;
plot(freq_sig/1e3,10*log10(abs(sig_fft(2,fidx)).^2)-sig_fft_max(2),...
     'k--','linewidth',3);
hold on
plot(freq_sig/1e3,10*log10(abs(sig_fft(3,fidx)).^2)-sig_fft_max(3),...
     'k-','linewidth',1);
axis([20 80 -40 5]);
legend('Ideal tx','Actual tx','location','south');
set(gca,'xtick',20:20:80,'ytick',-40:20:0);
set(gca,'fontsize',12,'ticklength',[1 1]*0.02);
xlabel('Frequency (kHz)','fontsize',14);
ylabel('Arbitrarily normalized unit (dB)','fontsize',14);
saveas(gcf,[SAVE_DIR,'/spectrum_ideal_actual.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/spectrum_ideal_actual.png'],'png');

figure;
plot(freq_sig/1e3,10*log10(abs(sig_fft(6,fidx)).^2)-sig_fft_max(2),...
     'k-','linewidth',1);
hold on
plot(freq_sig/1e3,10*log10(abs(sig_fft(7,fidx)).^2)-sig_fft_max(3),...
     'k--','linewidth',3);
axis([20 80 -40 5]);
legend('Narrowband tx','Broadband tx','location','south');
set(gca,'xtick',20:20:80,'ytick',-40:20:0);
set(gca,'fontsize',12,'ticklength',[1 1]*0.02);
xlabel('Frequency (kHz)','fontsize',14);
ylabel('Arbitrarily normalized unit (dB)','fontsize',14);
saveas(gcf,[SAVE_DIR,'/spectrum_nb_bb.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/spectrum_nb_bb.png'],'png');



% Autocorrelation function
figure;
plot(t_xcorr*1e3,sig_xcorr_env(5,:),'k','linewidth',1);
hold on
plot(t_xcorr*1e3,sig_xcorr_env(4,:),'k--','linewidth',3);
axis([-0.15 0.15 0 1.1]);
legend('Low Hann win','High Hann win','location','southeast');
set(gca,'xtick',-0.15:0.05:0.15,'ytick',0:0.5:1);
set(gca,'fontsize',12,'ticklength',[1 1]*0.02);
xlabel('Time (ms)','fontsize',14);
ylabel('Normalized envelope','fontsize',14);
saveas(gcf,[SAVE_DIR,'/acorr_hihann_lohann.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/acorr_hihann_lohann.png'],'png');

figure;
plot(t_xcorr*1e3,sig_xcorr_env(2,:),'k--','linewidth',3);
hold on
plot(t_xcorr*1e3,sig_xcorr_env(3,:),'k-','linewidth',1);
axis([-0.15 0.15 0 1.1]);
legend('Ideal tx','Actual tx','location','southeast');
set(gca,'xtick',-0.15:0.05:0.15,'ytick',0:0.5:1);
set(gca,'fontsize',12,'ticklength',[1 1]*0.02);
xlabel('Time (ms)','fontsize',14);
ylabel('Normalized envelope','fontsize',14);
saveas(gcf,[SAVE_DIR,'/acorr_ideal_actual.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/acorr_ideal_actual.png'],'png');

figure;
plot(t_xcorr*1e3,sig_xcorr_env(6,:),'k','linewidth',1);
hold on
plot(t_xcorr*1e3,sig_xcorr_env(7,:),'k--','linewidth',3);
axis([-0.15 0.15 0 1.1]);
legend('Narrowband tx','Broadband tx','location','southeast');
set(gca,'xtick',-0.15:0.05:0.15,'ytick',0:0.5:1);
set(gca,'fontsize',12,'ticklength',[1 1]*0.02);
xlabel('Time (ms)','fontsize',14);
ylabel('Normalized envelope','fontsize',14);
saveas(gcf,[SAVE_DIR,'/acorr_nb_bb.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/acorr_nb_bb.png'],'png');


