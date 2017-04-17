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
% 2012 05 21  use a fixed beampattern across the band
% 2012 08 13  use this code to do simluation WITHOUTH beampattern at all

function model_n_scat_no_bp(N,sampleN,gate_len,tx_opt,taper_opt,save_opt,save_dir,fname)
%INPUT
%   N           number of scatterers in the gate
%   sampleN     total number of realizations
%   gateLdist   model gate length
%   use_ideal   1-square chirp
%               2-ideal transmit signal (no system response)
%               3-actual transmit signal (with system response)
%   taper_opt   0-no taper
%               1-full Gaussian taper
%               2-HF Hann taper
%               3-LF Hann taper
%   save_opt    1-yes, 0-no
%   save_dir    folder to save the results too
%   fname       description of the simulation condition

%{
N = 5;
sampleN = 10;
gateLdist = 0.5;
save_opt=0;
taper_opt=0;  % taper transmit signal
use_ideal=0;  % 0-use ideal chirp for transmit signal
              % 1-use actual tx signal with system response
save_dir = '/mnt/storage/ECHO_STAT/20120521_bp_results';
%}


disp(['N=',num2str(N)]);
save_file = [fname,'_','N',num2str(N),...
             '_sampleN',num2str(sampleN),...
             '_gateLen',num2str(gate_len),'_noBP.mat'];
if taper_opt==1
    save_file = ['taper_',save_file];
end


%% Load in beampattern ir
BPIR = load('/mnt/storage/ECHO_STAT/20120223_bp_results/bpir_a0.054m_dtheta0.001pi_fmax10000kHz.mat');
dt = diff(BPIR.t_all(1:2));
bpirL = length(BPIR.t_all);
bpirHalfL = (length(BPIR.t_all)+1)/2;  % half length of the bpir

% Frame length parameters
c = 1500;  % sounds speed
gateL = round(gate_len/1500/dt);
%frameL = gateL+bpirL;
%frameL = gateL;
%t_frame = (0:frameL-1)*dt;

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
%[r,r_t] = xcorr(ones(1,frameL),yacorr);
%rL = length(r);

yHalfL = (length(yacorr)+1)/2;  % half length of the bpir
yL = length(yacorr);
rL = gateL+yL;

% Sample point
smpl_pt = yHalfL+(1:gateL)-1;
smpl_pt_mid = smpl_pt(round((length(smpl_pt)+1)/2));

% Calculation
s = zeros(1,sampleN);
%resp = zeros(sampleN,rL);
%resp_x = zeros(sampleN,length(smpl_pt));
%resp_x_env = zeros(sampleN,length(smpl_pt));
tic
parfor iS=1:sampleN

amp = raylrnd(1/sqrt(2)*ones(N,1));  % Rayleigh scatterer
time = round(rand(N,1)*gateL);
theta = rand_piston_angle(N)';  % angle in the beam
                                % theta = 10/180*pi;

resp_1smpl = zeros(1,rL);
for iN=1:N
    resp_temp = amp(iN)*...
        [zeros(1,time(iN)),yacorr,...
         zeros(1,rL-time(iN)-yL)];
    resp_1smpl = resp_1smpl + resp_temp;
    %figure;plot(resp_1smpl);
end

%resp(iS,:) = resp_1smpl;

%tmp = xcorr(resp_1smpl,yacorr);
%tmp = xcorr(yacorr,resp_1smpl);
%tmp2 = abs(hilbert(tmp));
tmp2 = abs(hilbert(resp_1smpl));
%resp_x(iS,:) = tmp(smpl_pt);
%resp_x_env(iS,:) = tmp2(smpl_pt);
s(iS) = tmp2(smpl_pt_mid);

end
disp('time to generate all samples')
toc


%resp_x = resp_x(:,smpl_pt);
%resp_x_env = resp_x_env(:,smpl_pt);

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
