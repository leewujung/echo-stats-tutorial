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


function model_n_scat_freqdep_bp_20120511(N,sampleN,gateLdist,use_ideal,taper_opt,save_opt,save_dir,fname)
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
N = 1;
sampleN = 10;
gateLdist = 0.5;
save_opt=0;
taper_opt=0;  % taper transmit signal
use_ideal=0;  % use ideal chirp for transmit signal
%}

disp(['N=',num2str(N)]);
save_file = [fname,'_','N',num2str(N),...
             '_sampleN',num2str(sampleN),...
             '_gateLdist',num2str(gateLdist),'_freqDepBP.mat'];
if taper_opt==1
    save_file = ['taper_',save_file];
end


% Load in beampattern ir
bp_file = '/mnt/storage/ECHO_STAT/20120223_bp_results/bpir_a0.054m_dtheta0.001pi_fmax10000kHz.mat';
BPIR = load(bp_file);
dt = diff(BPIR.t_all(1:2));
bpirL = length(BPIR.t_all);
bpirHalfL = (length(BPIR.t_all)+1)/2;  % half length of the bpir

% Frame length parameters
c = 1500;  % sounds speed
gateL = round(gateLdist/1500/dt);
frameL = gateL+bpirL;
t_frame = (0:frameL-1)*dt;

% Obtain replica with system response
if use_ideal==1   % Square chirp signal
    t_y = 0:dt:0.01;          % pulse length [sec]
    y = chirp(t_y,30e3,0.01,70e3);   % start/end freq [Hz]
elseif use_ideal==2  % Ideal transmit signal (no system resp)
    [y,t_y] = chirp_no_sys(120,100);
    y = y/max(y);
    % resample according to dt from bpir
    [p,q] = rat(diff(t_y(1:2))/dt);
    y = resample(y,p,q);
    t_y = ((1:length(y))-1)*dt;
elseif use_ideal==3  % Actual transmit signal (with system resp)
    [y,t_y] = chirp_w_sys(120,100);
    y = y/max(y);
    % resample according to dt from bpir
    [p,q] = rat(diff(t_y(1:2))/dt);
    y = resample(y,p,q);
    t_y = ((1:length(y))-1)*dt;
end

% windowing the chirp
if taper_opt==1  % full Gaussian taper
    winL = round(length(y)/2);
    win = gausswin(winL);
    win = [win(1:floor(winL/2))',ones(1,length(y)-winL),win(floor(winL/ ...
                                                      2)+1:end)'];
    y = y.*win;
elseif taper_opt==2  % HF Hann taper
    winL = round(length(y)/2);
    win = hann(winL);
    win = [zeros(1,length(y)/2),win(1:floor(winL/2))',ones(1,winL/2)];
    y = y.*win;
elseif taper_opt==3  % LF Hann taper
    winL = round(length(y)/2);
    win = hann(winL);
    win = [ones(1,winL/2),win(floor(winL/2)+1:end)',zeros(1,length(y)/2)];
    y = y.*win;
elseif taper_opt==4  % middle Hann taper
    winL = round(length(y)/4);
    win = hann(winL);
    win = [zeros(1,winL),win(1:winL/2)',ones(1,winL),win(winL/2+1:end)',zeros(1,winL)];
    y = y.*win;
elseif taper_opt==5  % upper Hann taper
    winL = round(length(y)/4);
    win = hann(winL);
    win = [zeros(1,winL*2),win(1:winL/2)',ones(1,winL),win(winL/2+1:end)'];
    y = y.*win;
elseif taper_opt==6  % lower Hann taper
    winL = round(length(y)/4);
    win = hann(winL);
    win = [win(1:winL/2)',ones(1,winL),win(winL/2+1:end)',zeros(1,winL*2)];
    y = y.*win;
elseif taper_opt==7  % middle Hann taper narrow1
    winL = round(length(y)/4);
    win = hann(winL);
    win = [zeros(1,winL*11/8),win(1:winL/2)',ones(1,winL/4),win(winL/2+1:end)',zeros(1,winL*11/8)];
    y = y.*win;
elseif taper_opt==77  % middle Hann taper narrow1 very narrow  1/4 y len
    winL = round(length(y)/4);
    win = hann(winL);
    win = [zeros(1,winL*3/2),win',zeros(1,winL*3/2)];
    y = y.*win;
elseif taper_opt==777  % middle Hann taper narrow1 very narrow  1/16 y len
    winL = round(length(y)/16);
    win = hann(winL);
    win = [zeros(1,winL*15/2),win',zeros(1,winL*15/2)];
    y = y.*win;
elseif taper_opt==8  % middle Hann taper narrow2 --> this one not used
    winL = round(length(y)/4);
    win = hann(winL);
    win = [zeros(1,winL*5/4),win(1:winL/2)',ones(1,winL/2),win(winL/2+1:end)',zeros(1,winL*5/4)];
    y = y.*win;
elseif taper_opt==9  % middle Hann taper wider
    winL = round(length(y)/4);
    win = hann(winL);
    win = [zeros(1,winL/2),win(1:winL/2)',ones(1,winL*2),win(winL/2+1:end)',zeros(1,winL/2)];
    y = y.*win;
elseif taper_opt==10  % middle Hann taper wider2
    winL = round(length(y)/4);
    win = hann(winL);
    win = [win(1:winL/2)',ones(1,winL*3),win(winL/2+1:end)'];
    y = y.*win;
elseif taper_opt==11  % upper Hann taper 20kHz band
    winL = round(length(y)/4);
    win = hann(winL);
    win = [zeros(1,winL),win(1:winL/2)',ones(1,winL*2),win(winL/2+1:end)'];
    y = y.*win;
elseif taper_opt==12  % lower Hann taper 20kHz band
    winL = round(length(y)/4);
    win = hann(winL);
    win = [win(1:winL/2)',ones(1,winL*2),win(winL/2+1:end)',zeros(1,winL)];
    y = y.*win;
end

[yacorr, t_yacorr] = xcorr(y);


% Get model result length
[r,r_t] = xcorr(ones(1,frameL),yacorr);
rL = length(r);

% Sample point
smpl_pt = (length(yacorr)+1)/2-bpirHalfL+(1:frameL+bpirL)-1;
smpl_pt_mid = smpl_pt(round((length(smpl_pt)+1)/2));

% Calculation
%resp = zeros(sampleN,frameL);
%resp_x = zeros(sampleN,length(smpl_pt));
%resp_x_env = zeros(sampleN,length(smpl_pt));

tic
parfor iS=1:sampleN
    
    %amp = raylrnd(1/sqrt(2)*ones(N,1));  % Rayleigh scatterer
    [amp,~,~] = prosph_3D_simulation(1/5,N,1);
    
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
        %figure;plot(resp_1smpl);
    end
    
    %resp(iS,:) = resp_1smpl;
    
    %tmp = xcorr(resp_1smpl,yacorr);
    tmp = xcorr(yacorr,resp_1smpl);
    tmp2 = abs(hilbert(tmp));
    %resp_x(iS,:) = tmp(smpl_pt);
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
param.gateLdist = gateLdist;
param.use_ideal = use_ideal;
param.taper_opt = taper_opt;
param.save_opt = save_opt;
param.save_dir = save_dir;
param.description = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n',...
                            'bp_file - file of which the beampattern impulse responses are stored',...
                            'N         - number of scatterer',...
                            'sampleN   - number of samples',...
                            'gateLdist - length of time gate',...
                           ['use_ideal - whether use ideal square chirp ',...
                              '(1-yes, 0-use actual tx with system response)'],...
                            'taper_opt - whether the tx is tapered or not (0-no, 1-yes)',...
                            'save_opt  - whether the results are saved (0-no, 1-yes)',...
                            'save_dir  - saved directory');


if save_opt==1
save([save_dir,'/',save_file],'N','sampleN','gateLdist',...
    's','smpl_pt','smpl_pt_mid','param');
%    's','smpl_pt','smpl_pt_mid');
end
%    'resp','resp_x','resp_x_env','s','smpl_pt','smpl_pt_mid');
