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
% 2012 11 10  extend this code to do mixed assemblages
% 2013 07 26  use the decimated transmit signal from the updated 'chirp_w_sys'


function model_n_scat_freqdep_bp_mixed(N,mix_r,sampleN,gate_len,tx_opt,taper_opt,save_opt,save_dir,fname,varargin)
%INPUT
%  N           number of scatterers in the gate
%              an array if mixed assemblage
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
%  mix_r       ratio between the components in mixed assemblages
%              =1 if simple aggregation, an array if mixed assemblage
%  varargin{1}  distribution for individual scatterer
%               0  - Rayliegh
%               >0 - aspect ratio for prolate spheroid
%               default=0 if not specified

%{
N = [10,10];
mix_r = [1,5];
sampleN = 10;
gate_len = 0.5;
save_opt = 0;
taper_opt = 0;  % taper transmit signal
tx_opt = 3;
%}

% Prep filename
if length(N)==1
    disp(['N=',num2str(N)]);
    save_file = [fname,'_','N',num2str(N),...
                 '_sampleN',num2str(sampleN),...
                 '_gateLen',num2str(gate_len),'_freqDepBP.mat'];
else
    nn = ['N_'];
    rr = ['r_'];
    for iN=1:length(N)
        nn = [nn,num2str(N(iN)),'_'];
        rr = [rr,num2str(mix_r(iN)),'_'];
    end
    disp([nn,rr]);
    save_file = [fname,'_',nn,rr,...
                 'sampleN',num2str(sampleN),...
                 '_gateLen',num2str(gate_len),'_freqDepBP.mat'];
end


%% Determine individual scatterer
if ~isempty(varargin)
    if varargin{1}==0 % Rayleigh scatterer
        indiv = 0;
        ar = [];
    else
        indiv = 1;
        ar = varargin{1};
    end
else
    indiv = 0;
    ar = [];
end


%% Load in beampattern ir
bp_folder = '/mnt/storage/broadband_code_current/bpir_pool';
bp_file = 'bpir_a0.054m_dtheta0.001pi_fmax10000kHz.mat';
BPIR = load([bp_folder,'/',bp_file]);
BPIR.dt = diff(BPIR.t_all(1:2));  % time spacing of bpir
bpirL = length(BPIR.t_all);
bpirHalfL = (length(BPIR.t_all)+1)/2;  % half length of the bpir

%% Obtain transmit signal/replica
[y,t_y] = gen_tx(tx_opt,BPIR.dt);
%[y_ori,t_y_ori] = gen_tx(tx_opt,BPIR.dt);
%y_ori_fft = fft(y_ori);
%freq_y_ori = 1/diff(t_y_ori(1:2))/(length(y_ori_fft)-1)*((1:(length(y_ori_fft)+1)/2)-1);

%% Resample bpir to fit the resolution of the replica
[p,q] = rat(BPIR.dt/diff(t_y(1:2)));
BPIR.bpir_all_rs = resample(BPIR.bpir_all,p,q);

%{
% resample according to dt from bpir
[p,q] = rat(diff(t_y(1:2))/dt);
y = resample(y,p,q);
t_y = ((1:length(y))-1)*dt;
y_fft = fft(y);
freq_y = 1/diff(t_y(1:2))/(length(y_fft)-1)*((1:(length(y_fft)+1)/2)-1);
%}

% window the chirp
win = win_chirp(taper_opt,y);
y = y.*win;
[yacorr, t_yacorr] = xcorr(y);

%% Frame length parameters
c = 1500;  % sounds speed
gateL = round(gate_len/1500/dt);
frameL = gateL+bpirL;
t_frame = (0:frameL-1)*dt;



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
resp_x = zeros(sampleN,length(smpl_pt)); % unmark this for resp_x
%resp_x_env = zeros(sampleN,length(smpl_pt));


tic
parfor iS=1:sampleN

    if indiv==0  % Rayleigh scatterer
        amp = [];
        for iN=1:length(N)
            amp = [amp;raylrnd(mix_r(iN)/sqrt(2)*ones(N(iN),1))];
        end
    else  % Randomly-rough prolate spheroid
        amp = [];
        for iN=1:length(N)
            [amp_tmp,~,~] = prosph_3D_simulation(ar,N(iN),1,mix_r(iN));
            amp = [amp,amp_tmp];
        end
    end
    
    time = round(rand(sum(N),1)*gateL + bpirL/2);
    
    theta = rand_piston_angle(sum(N))';  % angle in the beam
                                    % theta = 10/180*pi;
    [~,ind] = min(abs(repmat(theta,1,length(BPIR.theta))-...
                      repmat(BPIR.theta,sum(N),1)),[],2); % pick right bpir
    
    resp_1smpl = zeros(1,frameL);
    for iN=1:sum(N)
        resp_temp = amp(iN)*[zeros(1,time(iN)-bpirHalfL),...
                            BPIR.bpir_all(ind(iN),:),...
                            zeros(1,frameL-time(iN)-bpirHalfL+1)];
        resp_1smpl = resp_1smpl + resp_temp;
    end
    
    tmp = xcorr(yacorr,resp_1smpl);
    tmp2 = abs(hilbert(tmp));
    resp_x(iS,:) = tmp(smpl_pt);  % unmark for resp_x

    %resp_x_env(iS,:) = tmp2(smpl_pt);
    s(iS) = tmp2(smpl_pt_mid);
    
end
disp('time to generate all samples')
toc
%resp_x = resp_x(:,smpl_pt);
%resp_x_env = resp_x_env(:,smpl_pt);

param.bp_file = bp_file;
param.N = N;
param.mix_r = mix_r;
param.indiv = indiv;
param.ar = ar;
param.sampleN = sampleN;
param.gate_len = gate_len;
param.tx_opt = tx_opt;
param.taper_opt = taper_opt;
param.save_opt = save_opt;
param.save_dir = save_dir;


%{
if save_opt==1
    save([save_dir,'/',save_file],'N','mix_r','sampleN','gate_len',...
         's','smpl_pt','smpl_pt_mid','param');
end
%}


if save_opt==1
    save([save_dir,'/',save_file],'N','mix_r','sampleN','gate_len',...
         'resp_x','s','smpl_pt','smpl_pt_mid','param');
end
