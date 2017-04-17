% 2013 07 27  Test if time and freq results are interchangeable
%             this is done to see if can speed up the broadband code by
%             doing things in the freq domain only. If yes, this will
%             also make incorporating fish scattering response easier

tx_opt = 3;

theta = 1.35/180*pi;
fmax = 1.5e6;
a = 0.054;
f_bp = 0:500:1.5e6;
bpf = bpf_2way_fcn(theta,f_bp,a);
[t_bp,bpt] = bpir_2way_fcn(theta,fmax,a);
dt = diff(t_bp(1:2));

[y,t_y] = gen_tx(tx_opt);
y_fft = fft(y);
freq_y = 1/diff(t_y(1:2))/(length(y_fft)-1)*((1:(length(y_fft)+1)/2)-1);

y_fft_acorr = conj(y_fft).*y_fft;

y_acorr = xcorr(y);
t_y_acorr = (0:length(y_acorr)-1)*diff(t_y(1:2));
y_acorr_fft = fft(y_acorr);
freq_y_acorr = 1/diff(t_y(1:2))/(length(y_acorr_fft)-1)*((1:(length(y_acorr_fft)+1)/2)-1);

%% Beampattern modification in the freq domain
bpf_y = interp1(f_bp,bpf,freq_y);
e_fft = y_fft_acorr(1:length(freq_y)).*bpf_y;
e_fft_ifft = ifftshift(ifft([e_fft,fliplr(conj(e_fft(2:end)))]));
e_fft_ifft_env = abs(hilbert(e_fft_ifft));
t_e_fft_ifft = (0:length(e_fft_ifft)-1)/(2*max(freq_y));
[~,idx_e_fft_ifft] = max(e_fft_ifft_env);

%% Beampattern modification in the time domain
% resample the beampattern IR to match y_acorr
[p,q] = rat(dt/diff(t_y(1:2)));
bpt_y = resample(bpt,p,q);
t_bpt_y = (0:length(bpt_y)-1)*dt*q/p;
et = conv(y_acorr,bpt_y);
et_fft = fft(et);
et_env = abs(hilbert(et));
t_et = (0:length(et)-1)*dt*q/p;
[~,idx_et] = max(et_env);
freq_et = 1/diff(t_y(1:2))/(length(et_fft)-1)*((1:(length(et_fft)+1)/2)-1);

% resample y_acorr to match the beampattern IR
[p,q] = rat(diff(t_y(1:2))/dt);
y_acorr_rs = resample(y_acorr,p,q);
t_y_acorr_rs = (0:length(y_acorr_rs)-1)*dt;
et2 = conv(y_acorr_rs,bpt);
et2_fft = fft(et2);
et2_env = abs(hilbert(et2));
[~,idx_et2] = max(et2_env);
freq_et2 = 1/dt/(length(et2_fft)-1)*((1:(length(et2_fft)+1)/2)-1);
t_et2 = (0:length(et2)-1)*dt;

%% Plot
% Freq domain check
figure;
plot(freq_y,20*log10(abs(e_fft(1:length(freq_y)))));
hold on
plot(freq_et,20*log10(abs(et_fft(1:length(freq_et)))),'r');
plot(freq_et2,20*log10(abs(et2_fft(1:length(freq_et2)))),'g');
legend({'freq domain op',...
        'time domain op: resample bt to match y',...
        'time fomain op: resample y to match bt'});
xlim([0e3 150e3]);

figure;
plot(t_e_fft_ifft-t_e_fft_ifft(idx_e_fft_ifft),...
     e_fft_ifft_env/max(e_fft_ifft_env));
hold on
plot(t_et-t_et(idx_et),et_env/max(et_env),'r');
plot(t_et2-t_et2(idx_et2),et2_env/max(et2_env),'g');
legend({'freq domain op',...
        'time domain op: resample bt to match y',...
        'time fomain op: resample y to match bt'});
xlim([-2e-4 2e-4]);
