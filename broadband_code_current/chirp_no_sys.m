function [y_deci,t_y_deci] = chirp_no_sys(towdepth,towpower)
% Return the transmitted chirp without system response
% will use it as the MF replica
% 2013 07 26  Decimate the transmit signal by 10x

% Ideal TX/RX signal
%towdepth = 120;
%towpower = 100;
cal_folder = '/mnt/storage/calibration_2008/Cal_mat_files';
cal_file = [num2str(towdepth),'m_AMLO_',num2str(towpower),'.mat'];
load([cal_folder,'/',cal_file],'data');
t_y = data.tx0.tx;
y = data.tx0.yxmit;  % transmit signal, very fine spacing
y_fft = fft(y);
freq_y = 1/diff(t_y(1:2))/(length(y_fft)-1)*((1:(length(y_fft)+1)/2)-1);

% Decimate the actual TX signal
deci_r = 9;
y_deci = decimate(y,deci_r);
y_deci_fft = fft(y_deci);
t_y_deci = (0:length(y_deci)-1)*diff(t_y(1:2))*deci_r;
freq_y_deci = 1/diff(t_y_deci(1:2))/(length(y_deci_fft)-1)*((1:(length(y_deci_fft)+1)/2)-1);
