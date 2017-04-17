% 2013 07 30  Test for generating the time domain impulse after modifying
%             the transmit signal with the fish scattering response
% 2013 08 01  Test for doing Hilbert transform on the time domain fish
%             scattering response --> see if there's imaginary part in
%             the raw response that produces warnings in Hilbert transform

tx_opt = 3;

[y,t_y] = gen_tx(tx_opt);
y_fft = fft(y);
freq_y = 1/diff(t_y(1:2))/(length(y_fft)-1)*((1:(length(y_fft)+1)/2)-1);
dt = 1/(2*freq_y(end));

Rss = conj(y_fft).*y_fft;
Rss = Rss(1:length(freq_y)).';

%% Beampattern modification in the freq domain
theta = (0:1:90)/180*pi;
fmax = 1.5e6;
a = 0.054;
f_bp = 0:500:1.5e6;
bpf = bpf_2way_fcn(theta,f_bp,a);
H_bp = interp1(f_bp,bpf,freq_y)';

record = zeros(size(theta));
for iT=1:length(theta)
%% Fish scattering response modification in the freq domain
fish_folder = '/mnt/storage/broadband_code_current/fish_info';
fish_file = 'fish_scat_response_angle-90to90deg_len19to29cm.mat';
FISH = load([fish_folder,'/',fish_file]);
FISH.fbs_y = interp1(FISH.freq_fish,FISH.fbs_len_angle,freq_y);

H_scat = repmat(Rss,1,size(FISH.fbs_y,2)).*FISH.fbs_y.*repmat(H_bp(:,iT),1,size(FISH.fbs_y,2));
H_scat(isnan(H_scat)) = 0;

h_scat = ifftshift(ifft([H_scat;flipud(conj(H_scat(2:end,:)))]));
if ~isreal(h_scat)
    record(iT) = 1;
end
h_scat_env = abs(hilbert(h_scat));

t_h = (0:size(h_scat,1)-1)*dt;

end


%{
%% Plot
% Freq domain
figure;
plot(freq_y,20*log10(abs(Rss)),'k');
hold on
plot(freq_y,20*log10(abs(H_bp)),'r--');
plot(freq_y,20*log10(abs(H_fish)),'b--');
plot(freq_y,20*log10(abs(H_scat_bp_only)),'r');
plot(freq_y,20*log10(abs(H_scat_fish_only)),'b');
plot(freq_y,20*log10(abs(H_scat)),'c--');
legend('Rss','bp','fish',...
       'scat bp only','scat fish only','scat all');

% Time domain
figure;
plot(t_h,h_scat,'c');
hold on
plot(t_h,h_scat_bp_only,'r');
plot(t_h,h_scat_fish_only,'b');
legend('scat all','scat bp only','scat fish only');
%}



% Fish freq response overview
figure;
imagesc(freq_y/1e3,FISH.angle/pi*180,...
        20*log10(abs(FISH.fbs_y(:,1:181)))');
xlabel('Frequency (kHz)');
ylabel('Angle from normal incidence (deg)');
caxis([-80 -30]);
colorbar
title('Fish frequency response');


