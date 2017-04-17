% Compute the xcorr result of the new bpir and the chirp

load('/mnt/storage/ECHO_STAT/20111226_bpir_pool/bpir_dtheta_0.01pi.mat')
%bpir_all = [zeros(100,10000),bpir_all,zeros(100,10000)];

dt = mean(diff(t_all));
t_y = 0:dt:0.001;          % pulse length [sec]
y = chirp(t_y,30e3,0.01,70e3);   % start/end freq [Hz]
[yacorr, t_yacorr] = xcorr(y);
%yacorr = [zeros(1,50000),yacorr,zeros(1,50000)];

[tmp,t_xcorr] = xcorr(yacorr,bpir_all(end,:),'coeff');

bp_xcorr = zeros(size(bpir_all,1),length(tmp));
bp_xcorr_env = zeros(size(bpir_all,1),length(tmp));
for iT=1:size(bpir_all,1)
    bp_xcorr(iT,:) = xcorr(yacorr,bpir_all(iT,:));
    bp_xcorr_env(iT,:) = abs(hilbert(bp_xcorr(iT,:)));
end