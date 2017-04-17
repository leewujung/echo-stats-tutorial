% 2014 06 25  Transmit signal characteristics: spectrum,
%             beampattern impulse response, and xcorr output


SAVE_DIR = '/mnt/storage/broadband_echopdf_ms_figs/fig_tx_bpir';
if ~exist(SAVE_DIR,'dir')
    mkdir(SAVE_DIR);
end

%% Load signals
[y,t_y] = gen_tx(3);
y_fft = fft(y);
freq_y = 1/diff(t_y(1:2))/(length(y_fft)-1)*((1:(length(y_fft)+1)/2)-1);
yL = length(y);
yHalfL = round((yL+1)/2);
dt = 1/(2*freq_y(end));  % time step for the whole simulation
% autocorrelation
Rss = conj(y_fft).*y_fft;
Rss = Rss(1:length(freq_y)).';

%% Load beampattern impulse response
bpir_file = 'bpir_a0.054m_dtheta0.001pi_fmax10000kHz.mat';
BPIR = load([bp_folder,'/',bpir_file]);

%% Load beampattern transfer function
bp_folder = '/mnt/storage/broadband_code_current/bpir_bpf_pool';
bp_file = 'bpf_a0.054m_dtheta0.001pi_fmax1500kHz_df100Hz.mat';
BP = load([bp_folder,'/',bp_file]);
BP.bp_y = interp1(BP.freq_bp,BP.bp,freq_y);

%% Convolve with beampattern transfer function
theta = [5,10,30]'/180*pi;  % angle in the beam
[~,ind] = min(abs(repmat(theta,1,length(BP.theta))-...
                  repmat(BP.theta,length(theta),1)),[],2); % pick the right bp
H_bp = BP.bp_y(:,ind);

H_scat = repmat(Rss,1,length(theta)).*H_bp;
h_scat = ifftshift(ifft([H_scat;flipud(conj(H_scat(2:end,:)))]),1);
h_scat_env = abs(hilbert(h_scat));

m = max(h_scat_env(:,1));
mididx = ceil(size(h_scat,1)/2);
idx = mididx+(-1000:1000);
t_h_scat = (-1000:1000)*dt;
[~,iii] = max(h_scat_env(:,1));

figure;
plot(freq_y/1e3,20*log10(abs(y_fft(1:length(freq_y))))-max(20*log10(abs(y_fft(1:length(freq_y))))));
ylim([-50 5]);
xlim([20 80]);
xlabel('Frequency (kHz)');
ylabel('Spectrum (dB)');
title('Transmit signal spectrum');
saveas(gcf,sprintf('%s/tx_spectrum.fig',SAVE_DIR),'fig');
print(sprintf('%s/tx_spectrum.png',SAVE_DIR),'-dpng');

figure;
plot(BPIR.t_all*1e3,BPIR.bpir_all(:,ind));
xlabel('Time (ms)');
ylabel('Arbitrary relative scale');
title('Beampattern impulse response');
saveas(gcf,sprintf('%s/bpir.fig',SAVE_DIR),'fig');
print(sprintf('%s/bpir.png',SAVE_DIR),'-dpng');

figure;
plot(t_h_scat*1e3-t_h_scat(iii),h_scat_env(idx,:)/m);
xlim([-0.2 0.2]);
ylim([0 1.1])
xlabel('Time (ms)');
ylabel('Normalized envelope');
title('Xcorr of Ryy and bpir');
saveas(gcf,sprintf('%s/bp_xcorr.fig',SAVE_DIR),'fig');
print(sprintf('%s/bp_xcorr.png',SAVE_DIR),'-dpng');


