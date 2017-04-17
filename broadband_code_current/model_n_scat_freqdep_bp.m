% 2011 11 22  simulate N Rayleigh-distributed scatterers
%             randomly located within a short window
%             this window corresponds to twice the width
%             of the widest beampattern ir
% 2012 02 14  make the frame length exact
%             sample at the middle of the frame
% 2012 02 16  test the effect of tapering (system response, etc.)
% 2012 02 23  update the bpir to the ones using the correct radius
% 2012 04 18  incorporte system response
% 2012 05 08  keep on workin to incorporate system response
% 2012 05 11  incorporate the use of the actual tx signal
% 2012 10 26  need to save the time series for noise addition
%             move tx signal generation and windowing to separate functions


function model_n_scat_freqdep_bp(N,sampleN,gate_len,tx_opt,taper_opt,save_opt,save_dir,fname,varargin)
%INPUT
%  N           number of scatterers in the gate
%  sampleN     total number of realizations
%  gate_len    model gate length
%  tx_opt      1-square chirp
%              2-ideal transmit signal (no system response)
%              3-actual transmit signal (with system response)
%  taper_opt   0-no taper
%              1-full Gaussian taper
%              2-HF Hann taper
%              3-LF Hann taper
%  save_opt    1-yes, 0-no
%  save_dir    folder to save the results too
%  fname       description of the simulation condition
%  varargin{1}  distribution for individual scatterer
%               1-Rayleigh
%               2-prolate spheroid
%  varargin{2}  if varargin{1}==2, varargin{2}=aspect ratio

%{
N = 10;
sampleN = 10;
gate_len = 0.5;
save_opt = 0;
taper_opt = 0;  % taper transmit signal
tx_opt = 3;
%}

disp(['N=',num2str(N)]);
save_file = [fname,'_','N',num2str(N),...
             '_sampleN',num2str(sampleN),...
             '_gateLen',num2str(gate_len),'_freqDepBP.mat'];
if taper_opt==1
    save_file = ['taper_',save_file];
end


%% Load in beampattern ir
bp_file = '/mnt/storage/ECHO_STAT/20120223_bp_results/bpir_a0.054m_dtheta0.001pi_fmax10000kHz.mat';
BPIR = load(bp_file);
dt = diff(BPIR.t_all(1:2));
bpirL = length(BPIR.t_all);
bpirHalfL = (length(BPIR.t_all)+1)/2;  % half length of the bpir


%% Frame length parameters
c = 1500;  % sounds speed
gateL = round(gate_len/1500/dt);
frameL = gateL+bpirL;
t_frame = (0:frameL-1)*dt;


%% Obtain replica with system response
[y,t_y] = gen_tx(tx_opt,dt);

% resample according to dt from bpir
[p,q] = rat(diff(t_y(1:2))/dt);
y = resample(y,p,q);
t_y = ((1:length(y))-1)*dt;

% window the chirp
win = win_chirp(taper_opt,y);
y = y.*win;
[yacorr, t_yacorr] = xcorr(y);


%% Get model result length
[r,r_t] = xcorr(ones(1,frameL),yacorr);
rL = length(r);


%% Sample point
smpl_pt = (length(yacorr)+1)/2-bpirHalfL+(1:frameL+bpirL)-1;
smpl_pt_mid = smpl_pt(round((length(smpl_pt)+1)/2));

mid = round((length(smpl_pt)+1)/2);
smpl_pt = smpl_pt(mid-2500:mid+2500);

% Calculation
%resp = zeros(sampleN,frameL);
%resp_x = zeros(sampleN,length(smpl_pt)); % unmark this for resp_x
%resp_x_env = zeros(sampleN,length(smpl_pt));

tic
parfor iS=1:sampleN
    
    if ~isempty(varargin)
        if varargin{1}==1   % Rayleigh scatterer
            amp = raylrnd(1/sqrt(2)*ones(N,1));
        elseif varargin{1}==2
            ar = varargin{2};  % randomly rough prolate spheroid
            [amp,~,~] = prosph_3D_simulation(ar,N,1);
        end
    else
        amp = raylrnd(1/sqrt(2)*ones(N,1)); % Rayleigh if not specified
    end

    time = round(rand(N,1)*gateL + bpirL/2);
    
    theta = rand_piston_angle(N)';  % angle in the beam
                                    % theta = 10/180*pi;
    [~,ind] = min(abs(repmat(theta,1,length(BPIR.theta))-...
                      repmat(BPIR.theta,N,1)),[],2); % pick right bpir
    
    resp_1smpl = zeros(1,frameL);
    for iN=1:N
        resp_temp = amp(iN)*[zeros(1,time(iN)-bpirHalfL),...
                            BPIR.bpir_all(ind(iN),:),...
                            zeros(1,frameL-time(iN)-bpirHalfL+1)];
        resp_1smpl = resp_1smpl + resp_temp;
    end
    
    tmp = xcorr(yacorr,resp_1smpl); % unmark below 3 lines for resp_x
    tmp2 = abs(hilbert(tmp));
    %resp_x(iS,:) = tmp(smpl_pt);  % unmark for resp_x

    %resp_x_env(iS,:) = tmp2(smpl_pt);
    s(iS) = tmp2(smpl_pt_mid);
    
end
disp('time to generate all samples')
toc
%resp_x = resp_x(:,smpl_pt);
%resp_x_env = resp_x_env(:,smpl_pt);

param.bp_file = bp_file;
param.N = N;
param.sampleN = sampleN;
param.gate_len = gate_len;
param.tx_opt = tx_opt;
param.taper_opt = taper_opt;
param.save_opt = save_opt;
param.save_dir = save_dir;


if save_opt==1
    save([save_dir,'/',save_file],'N','sampleN','gate_len',...
         's','smpl_pt','smpl_pt_mid','param');
end

%{
if save_opt==1
    save([save_dir,'/',save_file],'N','sampleN','gate_len',...
         'resp_x','s','smpl_pt','smpl_pt_mid','param');
end
%}