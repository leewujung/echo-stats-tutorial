function [y_new,t_y_new] = chirp_w_sys(towdepth,towpower)
% Calculate the transmitted chirp after system response
% will use it as the MF replica
% 2013 07 27  Originally use finely-spaced tx signal
%             now decimate the signal by 10x

% Actual TX signal
%towdepth = 120;
%towpower = 100;
cal_folder = '/mnt/storage/calibration_2008/Cal_mat_files';
cal_file = [num2str(towdepth),'m_AMLO_',num2str(towpower),'.mat'];
load([cal_folder,'/',cal_file],'data');
t_y = data.tx0.tx;
y = data.tx0.yxmit;  % transmit signal, very fine spacing
y_fft = fft(y);
freq_y = 1/diff(t_y(1:2))/(length(y_fft)-1)*((1:(length(y_fft)+1)/2)-1);
clear data

% Decimate the actual TX signal
deci_r = 9;
y_deci = decimate(y,deci_r);
y_deci_fft = fft(y_deci);
t_y_deci = (0:length(y_deci)-1)*diff(t_y(1:2))*deci_r;
freq_y_deci = 1/diff(t_y_deci(1:2))/(length(y_deci_fft)-1)*((1:(length(y_deci_fft)+1)/2)-1);

% Apply system response
A=load(['/mnt/storage/calibration_2008/cal_result_20120402/' ...
        'cal_result_depth_120_power_100.mat']);
freq_sys = A.result.freq_cal;
Hsys = A.result.Hsys_wL;   % system respone considering TL
clear A

idx = zeros(length(freq_sys),1);
idx(174:463) = 1;
idx = logical(idx);
Hsys(~idx) = nan;  % take out the band edge
%Hsys = Hsys-max(Hsys);  % old before 20130725
Hsys = Hsys-max(Hsys(idx));  % new at 20130725
Hsys = 10.^(Hsys/20);
Hsys(1) = 0;
Hsys(end) = 0;
idx2 = isnan(Hsys);
freq_sys(idx2) = [];
Hsys(idx2) = [];

Hsys_y = interp1(freq_sys,Hsys,freq_y_deci);

y_sys_fft = y_deci_fft(1:length(freq_y_deci)).*Hsys_y;
y_sys_fft(isnan(y_sys_fft))=0;
y_sys_fft = [y_sys_fft,fliplr(conj(y_sys_fft(2:end)))];
y_new = ifft(y_sys_fft);
t_y_new = ((1:length(y_new))-1)/(2*max(freq_y_deci));

%y_new_fft = fft(y_new);   % check if the same --> yes
