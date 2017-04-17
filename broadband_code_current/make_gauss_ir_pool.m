% 2012 02 27  make a pool of Gaussian windows

minCycleNum = 3;  % num of cycles in the shortest window
freq = 50e3;      % [Hz]
fs = 10e5;         % sampling freq [Hz]
dtheta = 0.001;
theta = (dtheta:dtheta:1)*pi/2;  % angle

maxCycleNum = minCycleNum*60;
width_fac = linspace(1,60,length(theta));


folder = '/mnt/storage/ECHO_STAT/20120229_amp_phase_theta_env';
fname = sprintf(['gauss_ir_dtheta%2.3fpi_freq%dkHz_fs%dkHz',...
                 '_minCycleNum%d_maxCycleNum%d.mat'],dtheta, ...
                freq/1e3,fs/1e3,minCycleNum,maxCycleNum);


%t1 = (1:round(minCycleNum/freq*fs))/fs;
%y1 = sin(2*pi*freq*t1);
tend = (1:round(maxCycleNum/freq*fs))/fs;
%yend = sin(2*pi*freq*tend);

% Want width = odd number
winL = length(tend);
if mod(winL,2)==0
    winL = length(tend)+1;
else
    winL = length(tend);
end
winHalfL = (winL-1)/2;

win_all = zeros(length(theta),winL);
y_all = zeros(length(theta),winL);
for iT=1:length(theta)
    t_pt = round(minCycleNum*width_fac(iT)/freq*fs);
    if mod(t_pt,2)==0
        t_tmp = -(t_pt/2):(t_pt/2);
    else
        t_tmp = -(t_pt-1)/2:(t_pt-1)/2;
    end
    sigma = length(t_tmp)/8;
    y = cos(2*pi*freq*t_tmp/fs);
    g = 1/sigma*exp(-t_tmp.^2/(2*sigma^2));
    if mod(t_pt,2)==0
        win_temp = [zeros(1,winHalfL-t_pt/2), g,...
                    zeros(1,winL-winHalfL-t_pt/2-1)];
        y_temp = [zeros(1,winHalfL-t_pt/2), g.*y,...
                    zeros(1,winL-winHalfL-t_pt/2-1)];
    else
        win_temp = [zeros(1,winHalfL-(t_pt-1)/2), g,...
                    zeros(1,winL-winHalfL-(t_pt-1)/2-1)];
        y_temp = [zeros(1,winHalfL-(t_pt-1)/2), g.*y,...
                    zeros(1,winL-winHalfL-(t_pt-1)/2-1)];
    end
    win_all(iT,:) = win_temp;
    y_all(iT,:) = y_temp;
end
t_all = (1:size(win_all,2))/fs;

save([folder,'/',fname],'t_all','win_all','theta','y_all',...
     'freq','fs','dtheta','minCycleNum','maxCycleNum','-MAT');
