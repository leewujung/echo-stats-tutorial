% 2011 12 07  perform random phasor sum for a particular freq sinusoid
%             the results are stored in '20111207_n_scat_narrowband'
% 2012 02 21  make it comparable to the broadband case
% 2014 08 31  revive to compare with bbechopdf

function nbechopdf(N,ns,gate_len,nb_freq,sdir,fname)
%N = 2;
%ns = 50000;
%gate_len = 0.5;
%nb_freq = 50e3;

% Parameters
c = 1500;    % sound speed [m]
dt = 5e-6;
a = 0.054;   % transducer diameter [m]


% Transmitting signal
yL = round(2*gate_len/c/dt); % > 2x gateLdist to ensure overlapping
yHalfL = floor((yL+1)/2);
t_y = (0:yL-1)*dt;
y = sin(2*pi*nb_freq*t_y);


% Gate length parameters
gateL = round(gate_len/c/dt);  % gate length in pt
frameL = gateL+yL;
frameHalfL = round((frameL+1)/2);
%frameHalfL = floor((frameL+1)/2);


% sample point
smpl_pt = (1:gateL)+yHalfL-1;
smpl_pt_mid = round((gateL+1)/2)+yHalfL-1;


%% Beampattern response
bp_folder = '/mnt/storage/broadband_code_current/bpir_bpf_pool';
bp_file = 'bpf_a0.054m_dtheta0.001pi_fmax1500kHz_df100Hz.mat';
BP = load([bp_folder,'/',bp_file]);
BP.bp_y = interp1(BP.freq_bp,BP.bp,nb_freq);


% Calculation
resp = zeros(ns,frameL);
resp_env = zeros(ns,length(smpl_pt));
s = zeros(ns,1);
resp_env2 = zeros(ns,length(smpl_pt));
s2 = zeros(ns,1);
parfor iS=1:ns
    amp = raylrnd(1/sqrt(2)*ones(N,1));  % Rayleigh scatterer

    theta = rand_piston_angle(N,[])';  % beampattern effect

    % direct calculation of bp
    fac = (2*pi*nb_freq/c)*a*sin(theta);
    bp = (2*besselj(1,fac)./(fac)).^2;

    % look up for bp from existing file
    %[~,ind] = min(abs(repmat(theta,1,length(BP.theta))-...
    %                  repmat(BP.theta,sum(N),1)),[],2); % pick the right bp
    %bp = BP.bp_y(ind);

    time = round(rand(N,1)*gateL + yHalfL);  % random time in pt

    resp_1smpl = zeros(1,frameL);
    for iN=1:N
        resp_temp = [zeros(1,time(iN)-yHalfL),...
                           amp(iN)*bp(iN)*y,...
                           zeros(1,frameL-time(iN)-yHalfL+1)];
        resp_1smpl = resp_1smpl + resp_temp;
    end
    
    tmp = abs(hilbert(resp_1smpl));
    s(iS) = tmp(smpl_pt_mid);
end


% Save modeling parameters
param.N = N;
param.ns = ns;
param.gate_len = gate_len;
param.nb_freq = nb_freq;
param.c = c;
param.dt = dt;
param.a = 0.054;


% Save file
if ~isempty(sdir)
    disp('saving file...');
    sfname = [fname,'_N',num2str(N),...
              '_sampleN',num2str(ns),...
              '_gateLen',num2str(gate_len),'_narrowband.mat'];
    save([sdir,'/',sfname],'param','s');
end

