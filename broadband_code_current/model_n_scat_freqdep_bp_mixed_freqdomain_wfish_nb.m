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
% 2013 08 02  further modification to make fish len distr adjustable
% 2013 08 06  make fish location in the beam adjustable
% 2013 08 07  try narrowband fish response
%             revise input
%                                                                     
function model_n_scat_freqdep_bp_mixed_freqdomain_wfish_nb(N,mix_r,ns,gate_len,tx,taper,sname,bpa,indiv,varargin)
%INPUT
%  N           number of scatterers in the gate
%              an array if mixed assemblage
%  ns     total number of realizations
%  gate_len    model gate length
%  tx      1-square chirp
%              2-ideal transmit signal (no system response)
%              3-actual transmit signal (with system response)
%  taper   0-no taper
%              1-full Gaussian taper
%              2-HF Hann taper
%              3-LF Hann taper
%  sopt    1-yes, 0-no
%  sdir    folder to save the results too
%  fname       description of the simulation condition
%  mix_r       ratio between the components in mixed assemblages
%              =1 if simple aggregation, an array if mixed assemblage
%  bpa     restricted angle in the beam
%              [] - entire half space
%              bpa - only within bpa [deg]
%  indiv    distribution for individual scatterer
%               0 - Rayliegh
%               1 - Fish
%               2 - prolate spheroid
%               3 - Narrowband fish response
%               default=0 if not specified
%  varargin{1}  if indiv=2: ar - aspect ratio of prolate spheroid 
%               if indiv=1 or 3: {len_bin,len_dist}, unit: m
%  varargin{2}  if indiv=1 or 3: {angle_mean,angle_std}, unit: deg
%  varargin{3}  if indiv=3: specified freq for narrowband fish response

%{
%N = [10,10];
%mix_r = [1,5];
N = [10];
mix_r = [1];
ns = 1e4;
gate_len = 0.5;
sopt = 1;
taper = 0;  % taper transmit signal
tx = 3;
sdir = '/mnt/storage/ECHO_STAT/20130802_test_pdfs';
indiv = 1;  % 0-Rayleigh, 1-fish
%fname = ['test_indiv',num2str(indiv)];
fname = ['test_indiv',num2str(indiv),'Fnearnormal'];
%fname = ['test_indiv',num2str(indiv),'22cmonly'];
%fname = ['test_indiv',num2str(indiv),'distr'];
%}


%% Prep filename
if length(N)==1
    disp(['N=',num2str(N)]);
    save_file = [fname,'_','N',num2str(N),...
                 '_sampleN',num2str(ns),...
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
                 'sampleN',num2str(ns),...
                 '_gateLen',num2str(gate_len),'_freqDepBP.mat'];
end

%% Transmit signal/replica
[y,t_y] = gen_tx(tx);

% window the chirp
win = win_chirp(taper,y);
y = y.*win;

y_fft = fft(y);
freq_y = 1/diff(t_y(1:2))/(length(y_fft)-1)*((1:(length(y_fft)+1)/2)-1);
yL = length(y);
yHalfL = round((yL+1)/2);
dt = 1/(2*freq_y(end));  % time step for the whole simulation
% autocorrelation
Rss = conj(y_fft).*y_fft;
Rss = Rss(1:length(freq_y)).';

%% Beampattern response
bp_folder = '/mnt/storage/broadband_code_current/bpir_bpf_pool';
bp_file = 'bpf_a0.054m_dtheta0.001pi_fmax1500kHz_df100Hz.mat';
BP = load([bp_folder,'/',bp_file]);
BP.bp_y = interp1(BP.freq_bp,BP.bp,freq_y);

%% Fish scattering info
% Load fish scattering model
fish_folder = '/mnt/storage/broadband_code_current/fish_info';
fish_file = 'fish_scat_response_angle-90to90deg_len19to29cm.mat';
FISH = load([fish_folder,'/',fish_file]);
FISH.fbs_y = interp1(FISH.freq_fish,FISH.fbs_len_angle,freq_y);
FISH.fbs_y(isnan(FISH.fbs_y)) = 0;

%% Determine individual scatterer
% initialize for parfor
ar = [];
len_bin = [];
len_dist = [];
angle_mean = [];
angle_std = [];
nb_freq = [];
% assign params
if indiv==0  % Rayleigh scatterer
                 %if indiv==0  % Rayleigh scatterer
    disp('Individual: Rayleigh scatterer');
elseif indiv==2  % prolate spheroid
    disp('Individual: prolate spheroid');
    if nargin==1
        ar = varargin{1};
    else
        ar = 5; % default aspect ratio =5
    end
elseif indiv==1 | indiv==3  % fish
    disp('Individual: fish');
    if nargin==12 | nargin==13
        if isempty(varargin{1})
            disp('Use default fish length distribution');
            fishL = load(['/mnt/storage/broadband_code_current/fish_info/fish_len_dist.mat']);
            len_bin = fishL.L_bin;
            len_dist = fishL.L_dist;
        else
            len_bin = varargin{1}(:,1);
            len_dist = varargin{1}(:,2);
        end
        if isempty(varargin{2})
            disp('Use default fish angle of orientation distribution');
            angle_mean = -13;
            angle_std = 10;
        else
            angle_mean = varargin{2}(1);
            angle_std = varargin{2}(2);
        end
    else
        disp('Use default fish length and angle of orientation distribution');
        fishL = load(['/mnt/storage/broadband_code_current/fish_info/fish_len_dist.mat']);
        len_bin = fishL.L_bin;
        len_dist = fishL.L_dist;
        angle_mean = -13;
        angle_std = 10;
    end
    if indiv==1
        disp('Broadband fish response');
    else
        if isempty(varargin{3})
            nb_freq = 50e3;
        else
            nb_freq = varargin{3};
        end
        disp(['Narrowband fish response: ',num2str(nb_freq/1e3),'kHz']); 
    end
end

% Adjust angles
angle_mean = angle_mean/180*pi;
angle_std = angle_std/180*pi;


%% Frame length parameters
c = 1500;  % sounds speed
gateL = round(gate_len/c/dt);
frameL = gateL+2*yL;
t_frame = (0:frameL-1)*dt;
mid_frame_pt = round((frameL+1)/2);

resp = zeros(frameL,ns);
%resp_env = zeros(frameL,ns);
s = zeros(1,ns);

%% Simulation loop
%H_fish = zeros(length(freq_y),sum(N));
tic
parfor iS=1:ns  % realization loop

    % Fish scattering response
    if indiv==0  % Rayleigh scatterer
        H_fish = [];
        for iN=1:length(N)
            H_fish = [H_fish, repmat(raylrnd(mix_r(iN)/sqrt(2),1,N(iN)),length(freq_y),1)];
            %H_fish = [H_fish, ones(length(freq_y),sum(N))];
        end
        %for iN=1:length(N)
        %    H_fish(:,cumsum([0 N(1:iN-1)])+(1:N(iN))) =...
        %        repmat(raylrnd(mix_r(iN)/sqrt(2),1,N(iN)),length(freq_y),1);
        %end

    elseif indiv==1 | indiv==3  % fish scattering response
        % randomly generate fish length
        [len,~] = discrete_rnd(len_bin,len_dist,sum(N));
        [~,len_idx] = min(abs(repmat(len,1,length(FISH.len))-...
                          repmat(FISH.len,sum(N),1)),[],2);
        % randomly generate the angle of orientation
        %angle = normrnd_truncated(angle_mean,angle_std,2,sum(N),[])'; % truncated normal distribution
        angle = unifrnd(angle_mean-angle_std,angle_mean+angle_std,sum(N),1); % uniform
        [~,angle_idx] = min(abs(repmat(angle,1,length(FISH.angle))-...
                          repmat(FISH.angle,sum(N),1)),[],2);
        % pick the right one from the pool
        idx = (len_idx-1)*length(FISH.angle)+angle_idx;
        if indiv==1
            H_fish = FISH.fbs_y(:,idx);
        else
            [~,nb_idx] = min(abs(freq_y-nb_freq));
            H_fish = repmat(FISH.fbs_y(nb_idx,idx),length(freq_y),1);
            H_fish(1,:) = 0;  % set freq=0 component=0
        end

    elseif indiv==2  % prolate spheroid scatterer
        H_fish = [];
        for iN=1:length(N)
            amp_tmp = prosph_3D_simulation(ar,N(iN),1,mix_r(iN));
            H_fish = [H_fish, repmat(amp_tmp,length(freq_y),1)];
        end
        %for iN=1:length(N)
        %    amp_tmp = prosph_3D_simulation(ar,N(iN),1,mix_r(iN));
        %    H_fish(:,cumsum([0 N(1:iN-1)])+(1:N(iN))) = repmat(amp_tmp,length(freq_y),1);
        %end
    end

    % Beampattern response
    theta = rand_piston_angle(sum(N),bpa)';  % angle in the beam
    [~,ind] = min(abs(repmat(theta,1,length(BP.theta))-...
                      repmat(BP.theta,sum(N),1)),[],2); % pick the right bp
    H_bp = BP.bp_y(:,ind);
    %H_bp = ones(length(freq_y),sum(N));

    % Assemble and ifft
    H_scat = repmat(Rss,1,sum(N)).*H_fish.*H_bp;
    h_scat = ifftshift(ifft([H_scat;flipud(conj(H_scat(2:end,:)))]),1);
    % Wrong use of ifftshift
    %h_scat = ifftshift(ifft([H_scat;flipud(conj(H_scat(2:end,:)))]));
    % Forgot to conjugate
    %h_scat = ifftshift(ifft([H_scat;flipud(H_scat(2:end,:))]));

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
    if ~isreal(resp_temp)
        disp('Error: Invalid time series with imaginary part');
    end
    resp_temp = abs(hilbert(resp_temp));
    %    resp_env(:,iS) = resp_temp;
    s(iS) = resp_temp(mid_frame_pt);

end % realization loop
disp('time to generate all samples')
toc

param.bp_file = bp_file;
param.fish_file = fish_file;
param.N = N;
param.mix_r = mix_r;
param.indiv = indiv;
param.sampleN = ns;
param.gate_len = gate_len;
param.tx = tx;
param.taper = taper;
param.sopt = sopt;
param.sdir = sdir;
if isempty(bpa)
    param.bp_location = 'entire half space';
else
    param.bp_location = ['restrict to within ',num2str(bpa),'deg polar angle'];
end
if indiv==1 | indiv==3
    param.fish_len_bin = len_bin;
    param.fish_len_dist = len_dist;
    param.fish_len_unit = 'm';
    param.fish_angle_mean = angle_mean/pi*180;
    param.fish_angle_std = angle_std/pi*180;
    param.fish_angle_unit = 'deg';
    if indiv==3
        param.fish_nb_freq = nb_freq;
    end
elseif indiv==2
    param.ar = ar;
end

if sopt==1
    save([sdir,'/',save_file],'param','s');
end

