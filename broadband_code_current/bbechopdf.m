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
% 2014 09 04  change the prolate spheroid part so that can
%             constrain the angle of orientation distribution
%             instead of always use all angles of orientation
%                                                                     
function bbechopdf(N,mix_r,ns,gate_len,tx,sdir,fname,bpa,indiv,varargin)
%INPUT
%1  N         number of scatterers in the gate
%            an array if mixed assemblage
%2  mix_r     ratio between the components in mixed assemblages
%            length(mix_r)=1 if simple aggregation, an array if mixed assemblage
%3  ns        total number of realizations
%4  gate_len  model gate length
%5  tx      [signal type, taper]
%          signal type: 1-square chirp, 2-ideal transmit signal (no
%          system response), 3-actual transmit signal (with system
%          response)
%          taper: 0-no taper, 1-full Gaussian taper, 2-HF Hann taper,
%          3-LF Hann taper
%6  sdir    folder to save the results to
%7  fname   description of the simulation condition
%8  bpa     restricted angle in the beam
%          [] - entire half space
%          bpa - only within bpa [deg]
%9  indiv   distribution for individual scatterer
%          0 - Rayliegh (default)
%          1 - Fish
%          2 - prolate spheroid
%          3 - Narrowband fish response
%10  varargin{1}  if indiv=2: [ar,len_bin,len_dist]
%                              ar-aspect ratio of prolate spheroid 
%                 if indiv=1 or 3: {len_bin,len_dist}, unit: m
%11  varargin{2}  if indiv=1-3: {angle_mean,angle_std}, unit: deg
%12  varargin{3}  if indiv=3: specified freq for narrowband fish response

%% Display N and mix_r
nn = ['N_'];
rr = ['r_'];
for iN=1:length(N)
    nn = [nn,num2str(N(iN)),'_'];
    rr = [rr,num2str(mix_r(iN)),'_'];
end
disp('-------------------');
disp([nn,rr]);


%% Transmit signal/replica
[y,t_y] = gen_tx(tx(1));
if length(tx)==2  % tapering
    win = win_chirp(taper,y);
    y = y.*win;
end

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


%% Fish scattering model
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
if nargin>12
    disp(['Invalid inputs: nargin=',num2str(nargin)]);
    return;
end

if nargin>9
    if indiv==0
        disp('Rayleigh scatterer');
        disp('Ignoring additional inputs');
    elseif indiv==2
        disp('prolate spheroid');
        if nargin>=10
            if ~isempty(varargin{1})
                if length(varargin{1})==1
                    ar = varargin{1};
                    disp('Use default fish length distribution');
                    fishL = load(['/mnt/storage/broadband_code_current/fish_info/fish_len_dist.mat']);
                    len_bin = fishL.L_bin;
                    len_dist = fishL.L_dist;
                elseif length(varargin{1})==3  % NEED TO FIX THIS
                                               % BUG SINCE len_bin
                                               % AND len_dist will
                                               % be vectors
                    ar       = varargin{1}(1);
                    len_bin  = varargin{1}(2);
                    len_dist = varargin{1}(3);
                    if size(len_bin,1)==1
                        len_bin = len_bin;
                        len_dist = len_dist;
                    end
                else
                    disp('Invalid input (prolate spheroid option)!');
                end
            else
                ar = 0.5;  % 2014/09/04 change from 5 to 0.5
                disp('Use default fish length distribution');
                fishL = load(['/mnt/storage/broadband_code_current/fish_info/fish_len_dist.mat']);
                len_bin = fishL.L_bin;
                len_dist = fishL.L_dist;
            end
            if ~isempty(varargin{2})
                if length(varargin{2})==1 & varargin{2}==-1  % use [0,2*pi] uniform distribution
                    angle_mean = NaN;
                    angle_std = NaN;
                else  % use prescribed distribution
                    angle_mean = varargin{2}(1);
                    angle_std  = varargin{2}(2);
                end
            else
                disp('Use default angle of orientation distribution');
                angle_mean = -13;
                angle_std = 10;
            end
            if nargin>11
                disp('Ignoring additional inputs');
            end
        end
    elseif indiv==1 | indiv==3
        disp('fish');
        if nargin<11
            disp('Need additional inputs!');
            return;
        else
            % assign fish length and angle distribution
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
            if indiv==1
                disp('broadband fish repsonse');
                if nargin>11
                    disp('Ignoring additional inputs');
                end
            end
            if indiv==3
                if nargin<12
                    disp('Need additional inputs!');
                    return;
                else
                    % assign nb_freq
                    if isempty(varargin{3})
                        nb_freq = 50e3;
                    else
                        nb_freq = varargin{3};
                    end
                    disp(['narrowband fish response at ',...
                          num2str(nb_freq/1e3),'kHz']);
                end
            end
        end
    else
        disp('Invalid individual type!');
        return;
    end
elseif nargin==9
    if indiv==0
        disp('Rayleigh scatterer');
    else
        disp('Need additional inputs!');
        return;
    end
else
    disp('Basic params not assigned!');
    return;
end


%% Adjust angles from [deg] to [rad]
angle_mean = angle_mean/180*pi;
angle_std = angle_std/180*pi;


%% Save parameters
param.bp_file = bp_file;
param.fish_file = fish_file;
param.N = N;
param.mix_r = mix_r;
param.indiv = indiv;
switch indiv
  case 0
    param.indiv_type = 'Rayleigh scatterer';
  case 1
    param.indiv_type = 'broadband fish scatterer';
  case 2
    param.indiv_type = 'prolate spheroid scatterer';
  case 3
    param.indiv_type = 'narrowband fish scatterer';
end
param.sampleN = ns;
param.gate_len = gate_len;
param.tx = tx;
if isempty(bpa)
    param.bp_location = 'entire half space';
else
    param.bp_location = ['restrict to within ',num2str(bpa),'deg polar angle'];
end
if indiv==1 | indiv==3
    param.fish_len_bin     = len_bin;
    param.fish_len_dist    = len_dist;
    param.fish_len_unit    = 'm';
    param.fish_angle_mean  = angle_mean/pi*180;
    param.fish_angle_std   = angle_std/pi*180;
    param.fish_angle_unit  = 'deg';
    if indiv==3
        param.fish_nb_freq = nb_freq;
    end
elseif indiv==2
    param.ar = ar;
    param.prosph_angle_mean = angle_mean/pi*180;
    param.prosph_angle_std  = angle_std/pi*180;
    param.prosph_len_bin    = len_bin;
    param.prosph_len_dist   = len_dist;
end
%param


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
        angle = normrnd_truncated(angle_mean,angle_std,2,sum(N),[])'; % truncated normal distribution
        %angle = unifrnd(angle_mean-angle_std,angle_mean+angle_std,sum(N),1); % uniform
        [~,angle_idx] = min(abs(repmat(angle,1,length(FISH.angle))-...
                          repmat(FISH.angle,sum(N),1)),[],2);
        % pick the right one from the pool
        idx = (len_idx-1)*length(FISH.angle)+angle_idx;
        if indiv==1   % broadband fish response
            H_fish = FISH.fbs_y(:,idx);
        else          % narrowband fish response
            [~,nb_idx] = min(abs(freq_y-nb_freq));
            H_fish = repmat(FISH.fbs_y(nb_idx,idx),length(freq_y),1);
            H_fish(1,:) = 0;  % set freq=0 component=0
        end

    elseif indiv==2  % prolate spheroid scatterer
        H_fish = [];
        for iN=1:length(N)
            %amp_tmp = prosph_3D_simulation(ar,N(iN),1,mix_r(iN));  % all angles of orientation
            amp_tmp = prosph_amp(ar,N(iN),1,mix_r(iN),...
                                 [angle_mean,angle_std],[len_bin',len_dist']);
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


%% Save file
if ~isempty(sdir)
    disp('saving file...');
    sfname = [fname,'_',nn,rr,...
              'sampleN',num2str(ns),...
              '_gateLen',num2str(gate_len),'_freqDepBP.mat'];
    save([sdir,'/',sfname],'param','s');
    %save([sdir,'/',sfname],'param');
end
