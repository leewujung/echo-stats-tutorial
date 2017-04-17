% 2011 12 07  perform random phasor sum for a particular freq sinusoid
%             the results are stored in '20111207_n_scat_narrowband'
% 2012 02 21  make it comparable to the broadband case


function narrowband_n_scat(N,sampleN,gateLdist,taper_opt,save_opt)
%N = 2;
%sampleN = 1e4;
%gateLdist = 0.5;  % gate length [m]
%taper_opt = 0;
%save_opt = 1;

% Parameters
freq = 50e3;  % [Hz]
c = 1500;    % sound speed [m]
dt = 5e-8;
a = 0.054;   % transducer diameter [m]

% Dir
disp(['N=',num2str(N)]);
save_dir = '/mnt/storage/ECHO_STAT/20120106_bp_results';
if taper_opt==1
    save_file = ['NB_freq_',num2str(freq/1e3),'kHz_taper_N',num2str(N),...
                 '_sampleN_',num2str(sampleN),...
                 '_gateLdist_',num2str(gateLdist),'.mat'];
else
    save_file = ['NB_freq_',num2str(freq/1e3),'kHz_N',num2str(N),...
                 '_sampleN_',num2str(sampleN),...
                 '_gateLdist_',num2str(gateLdist),'.mat'];
end

% Transmitting signal
yL = round(2*gateLdist/c/dt); % > 2x gateLdist to ensure overlapping
yHalfL = floor((yL+1)/2);
t_y = (0:yL-1)*dt;
y = sin(2*pi*freq*t_y);

% Gate length parameters
gateL = round(gateLdist/c/dt);  % gate length in pt
frameL = gateL+yL;
frameHalfL = floor((frameL+1)/2);

% sample point
smpl_pt = (1:gateL)+yHalfL-1;
smpl_pt_mid = round((gateL+1)/2)+yHalfL-1;

% Calculation
resp = zeros(sampleN,frameL);
resp_env = zeros(sampleN,length(smpl_pt));
s = zeros(sampleN,1);
parfor iS=1:sampleN
    disp(['iS=',num2str(iS)]);

    amp = raylrnd(1/sqrt(2)*ones(N,1));  % Rayleigh scatterer

    theta = rand_piston_angle(N)';  % beampattern effect
    fac = (2*pi*freq/c)*a*sin(theta);
    bp = (2*besselj(1,fac)./(fac)).^2;

    time = round(rand(N,1)*gateL + yHalfL);  % random time in pt

    resp_1smpl = zeros(1,frameL);
    for iN=1:N
        resp_temp = [zeros(1,time(iN)-yHalfL),...
                           amp(iN)*bp(iN)*y,...
                           zeros(1,frameL-time(iN)-yHalfL+1)];
        resp_1smpl = resp_1smpl + resp_temp;
    end
    
    resp(iS,:) = resp_1smpl;
    tmp = abs(hilbert(resp_1smpl));
    resp_env(iS,:) = tmp(smpl_pt);
    s(iS) = tmp(smpl_pt_mid);
end

if save_opt==1
save([save_dir,'/',save_file],'N','sampleN','gateLdist',...
    'resp','s','smpl_pt','smpl_pt_mid');
end
