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
% 2013 07 27  use the decimated transmit signal from the updated 'chirp_w_sys'
% 2013 07 29  do all calculation in the freq domain so that it's easier
%             to incorporate both the beampattern and fish scattering response


function model_n_scat_freqdep_bp_mixed_freqdomain(N,mix_r,sampleN,gate_len,tx_opt,taper_opt,save_opt,save_dir,fname,varargin)
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
%N = [10,10];
%mix_r = [1,5];
N = [1];
mix_r = [1];
sampleN = 5e4;
gate_len = 0.5;
save_opt = 0;
taper_opt = 0;  % taper transmit signal
tx_opt = 3;
fname = 'test';
save_dir = 'test';
%}

disp('Old routine');

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
%indiv = 0;
%ar = [];
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

%% Obtain transmit signal/replica
[y,t_y] = gen_tx(tx_opt);

% window the chirp
win = win_chirp(taper_opt,y);
y = y.*win;

y_fft = fft(y);
freq_y = 1/diff(t_y(1:2))/(length(y_fft)-1)*((1:(length(y_fft)+1)/2)-1);
yL = length(y);
yHalfL = round((yL+1)/2);
dt = 1/(2*freq_y(end));  % time step for the whole simulation
% autocorrelation
Rss = conj(y_fft).*y_fft;
Rss = Rss(1:length(freq_y));

%% Load in beampattern response
bp_folder = '/mnt/storage/broadband_code_current/bpir_bpf_pool';
bp_file = 'bpf_a0.054m_dtheta0.001pi_fmax1500kHz_df100Hz.mat';
BP = load([bp_folder,'/',bp_file]);
BP.bp_y = interp1(BP.freq_bp,BP.bp,freq_y);


%% Frame length parameters
c = 1500;  % sounds speed
gateL = round(gate_len/c/dt);
frameL = gateL+2*yL;
t_frame = (0:frameL-1)*dt;
mid_frame_pt = round((frameL+1)/2);

resp = zeros(frameL,sampleN);
%resp_env = zeros(frameL,sampleN);
s = zeros(1,sampleN);

tic
parfor iS=1:sampleN

    % Fish scattering response
    if indiv==0  % Rayleigh scatterer
        H_fish = [];
        for iN=1:length(N)
            H_fish = [H_fish, repmat(raylrnd(mix_r(iN)/sqrt(2),1,N(iN)),length(freq_y),1)];
        end
    else  % fish scattering response
    end

    % Beampattern response
    theta = rand_piston_angle(sum(N),[])';  % angle in the beam
                                    % theta = 10/180*pi;
    [~,ind] = min(abs(repmat(theta,1,length(BP.theta))-...
                      repmat(BP.theta,sum(N),1)),[],2); % pick the right bp
    H_bp = BP.bp_y(:,ind);

    % Assemble and ifft
    H_scat = repmat(Rss.',1,sum(N)).*H_fish.*H_bp;
    h_scat = ifftshift(ifft([H_scat;flipud(conj(H_scat(2:end,:)))]),1);

    % Delay (location of fish)
    % need to do in the time domain since phase variation > 2*pi
    time = round(rand(sum(N),1)*gateL + yL);
    
    % Time domain impulse summation
    resp_temp = zeros(frameL,sum(N));
    for iN=1:sum(N)
        resp_temp(:,iN) = [zeros(time(iN)-yHalfL,1);h_scat(:,iN);...
                           zeros(frameL-time(iN)-yHalfL+1,1)];
    end
    resp_temp = sum(resp_temp,2);
    resp(:,iS) = resp_temp;
    resp_temp = abs(hilbert(resp_temp));
    %    resp_env(:,iS) = resp_temp;
    s(iS) = resp_temp(mid_frame_pt);

end
disp('time to generate all samples')
toc

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

if save_opt==1
    save([save_dir,'/',save_file],'param','s');
end

