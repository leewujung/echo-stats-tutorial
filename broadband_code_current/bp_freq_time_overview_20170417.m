% 2013 07 30  Test for generating the time domain impulse after modifying
%             the transmit signal with the fish scattering response
% 2013 08 01  Test for doing Hilbert transform on the time domain fish
%             scattering response --> see if there's imaginary part in
%             the raw response that produces warnings in Hilbert transform
% 2017 04 17  Beampattern response in frequency domain for 3 deg beam
%             use ka=44.2511 from 'fig_12_pb_ka_ka_num.mat'

addpath '~/Dropbox/0_CODE'/MATLAB/saveSameSize/

base_path = '/Volumes/wjlee_apl_2/echo_stat_tutorial/echo_stat_figs/';

% Make save path
str = strsplit(mfilename('fullpath'),'/');
str = str{end};
save_path = fullfile(base_path,str);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Transmit signal
% One octave band centered at 50 kHz
fctr = 50e3;
fd = 2^(1/2);
fup = fctr*fd;
flo = fctr/fd;

% Generate signal
t_y = 0:2.5e-6:0.01;    % pulse length [sec]
y = chirp(t_y,flo,0.01,fup);   % start/end freq [Hz]

y_fft = fft(y);
freq_y = 1/diff(t_y(1:2))/(length(y_fft)-1)*((1:(length(y_fft)+1)/2)-1);
dt = 1/(2*freq_y(end));

% Beampattern response
bp_folder = fullfile(base_path,'make_bpf_pool_20170417');
bp_file = 'bpf_a0.211m_dtheta0.010pi_fmax1500kHz_df100Hz.mat';
BP = load([bp_folder,'/',bp_file]);
BP.bp_y = interp1(BP.freq_bp,BP.bp,freq_y);

% Cross-correlation
Rss = conj(y_fft).*y_fft;
Rss = Rss(1:length(freq_y)).';

H_scat = repmat(Rss,1,size(BP.bp_y,2)).*BP.bp_y;
H_scat(isnan(H_scat)) = 0;

h_scat = ifftshift(ifft([H_scat;flipud(conj(H_scat(2:end,:)))]),1);
h_scat_env = abs(hilbert(h_scat));

t_h = (0:size(h_scat,1)-1)*dt;


% beampattern freq response overview
figure;
imagesc(BP.theta/pi*180,freq_y/1e3,...
        20*log10(abs(BP.bp_y)));
ylabel('Frequency (kHz)');
xlabel('Polar angle (deg)');
colorbar
title('Beampattern frequency response');
saveas(gcf,fullfile(save_path,'bp_freq_resp_overview.fig'),'fig');
saveSameSize_150(gcf,'file',fullfile(save_path,'bp_freq_resp_overview.png'),...
    'format','png');

% Fish time domain response overview
figure;
imagesc(t_h*1e3,BP.theta/pi*180,h_scat_env');
xlabel('Time (ms)');
ylabel('Polar angle (deg)');
colorbar
title('Beampattern time domain response (with Rss), linear color');
xlim([5.1 5.6]);
saveas(gcf,fullfile(save_path,'bp_time_resp_overview_linear.fig'),'fig');
saveSameSize_150(gcf,'file',fullfile(save_path,'bp_time_resp_overview_linear.png'),...
    'format','png');

figure;
imagesc(t_h*1e3,BP.theta/pi*180,20*log10(h_scat_env'));
xlabel('Time (ms)');
ylabel('Polar angle (deg)');
colorbar
c = caxis;
caxis([-60 c(2)])
title('Beampattern time domain response (with Rss), log color');
xlim([5.1 5.6]);
saveas(gcf,fullfile(save_path,'bp_time_resp_overview_log.fig'),'fig');
saveSameSize_150(gcf,'file',fullfile(save_path,'bp_time_resp_overview_log.png'),...
    'format','png');

