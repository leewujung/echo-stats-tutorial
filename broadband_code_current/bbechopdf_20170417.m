function [s,param] = bbechopdf_20170417(N,mix_r,num_sample,param)
% INPUT
%   N                       number of scatterers in the gate
%                           an array if mixed assemblage
%   mix_r                   ratio between the components in mixed assemblages
%                           length(mix_r)=1 if simple aggregation, an array if mixed assemblage
%   num_sample              total number of realizations
%   gate_len                model gate length
%   param.tx_opt            1-square chirp
%                           2-ideal transmit signal (no system response)
%                           3-actual transmit signal (with system response)
%   param.tx_opt_taper      0-no taper
%                           1-full Gaussian taper
%                           2-HF Hann taper
%                           3-LF Hann taper
%   param.bpa               restricted angle in the beam
%                           [] - entire half space
%                           bpa - only within bpa [deg]
%   param.nbwb              use narrowband (nb) or broadband (wb) beampatern and scatterer response
%   param.scatterer.type    distribution for individual scatterer
%                           choices: rayl, point, rayl, prosph, fish
%   param.scatterer.ar
%   param.scatterer.len_bin
%   param.scatterer.len_dist
%   param.scatterer.len_unit
%   param.scatterer.angle_mean
%   param.scatterer.angle_std
%   param.scatterer.angle_unit 
%
% Wu-Jung Lee | leewujung@gmail.com
% 
% 2017 04 17  Overhaul old bbechopdf function to use struct as input
%
% STILL NEED WORK:
%   -- prolate spheroid option
%   -- get_tx function, now generate an octave bandwidth chirp centered at
%      50 kHz
%   -- restrained beampattern angle: variable 'bpa'


%% Display N and mix_r
nn = 'N_';
rr = 'r_';
for iN=1:length(N)
    nn = [nn,num2str(N(iN)),'_'];
    rr = [rr,num2str(mix_r(iN)),'_'];
end
disp('-------------------');
disp([nn,rr]);


%% Save params from input
param.N = N;
param.mix_r = mix_r;
param.num_sample = num_sample;



%% Determine individual scatterer
% Assign default values
if strcmp(param.scatterer.type,'prosph')
    if ~isfield(param.scatterer,'ar')
        param.scatterer.ar = 0.5;
    end
    if ~isfield(param.scatterer,'sph_rot_opt')
        param.scatterer.sph_rot_opt = '3D';
    end
end

if strcmp(param.scatterer.type,'fish')
    if ~isfield(param.scatterer,'len_bin')
        fishL = load(fullfile(param.fish_len_path,param.fish_len_file));
        param.scatterer.len_bin = fishL.L_bin;
        param.scatterer.len_dist = fishL.L_dist;
        param.scatterer.len_unit = 'm';
    end
    if ~isfield(param.scatterer,'angle_mean')
        param.scatterer.angle_mean = -13;
        param.scatterer.angle_std = 10;
        param.scatterer.angle_unit  = 'deg';
    end
end

if ~isfield(param,'bpa')
    param.bpa = [];
end

if ~isfield(param,'nb_freq')
    param.nb_freq = 50e3;  % dummy for parfor even if in 'wb' mode
end



%% Transmit signal/replica
% One octave band centered at 50 kHz
fctr = 50e3;
fd = 2^(1/2);
fup = fctr*fd;
flo = fctr/fd;

% Generate signal
if strcmp(param.nbwb,'wb')
    t_y = 0:2.5e-6:0.01;    % pulse length [sec]
    y = chirp(t_y,flo,0.01,fup);   % start/end freq [Hz]    
    y_fft = fft(y);
    freq_y = 1/diff(t_y(1:2))/(length(y_fft)-1)*((1:(length(y_fft)+1)/2)-1);
    dt = 1/(2*freq_y(end));
    yL = length(y);
    % autocorrelation
    Rss = conj(y_fft).*y_fft;
    Rss = Rss(1:length(freq_y)).';
%     [~,nb_idx] = min(abs(freq_y-param.nb_freq));

elseif strcmp(param.nbwb,'nb')
    dt = 2.5e-6;
    freq_y = param.nb_freq;  % dummy for parfor
    yL = round(2*param.gate_len/param.c/dt); % > 2x gateLdist to ensure overlapping
    t_y = (0:yL-1)*dt;
    y = sin(2*pi*param.nb_freq*t_y);
    Rss= [];  % dummy for parfor
end

yHalfL = floor((yL+1)/2);

% [y,t_y] = gen_tx(param.tx_opt);
% if isfield(param,'tx_taper')  % tapering
%     win = win_chirp(param.tx_taper,y);
%     y = y.*win;
% end

% y_fft = fft(y);
% freq_y = 1/diff(t_y(1:2))/(length(y_fft)-1)*((1:(length(y_fft)+1)/2)-1);
% yL = length(y);
% % yHalfL = round((yL+1)/2);
% yHalfL = floor((yL+1)/2);
% dt = 1/(2*freq_y(end));  % time step for the whole simulation





%% Beampattern response
BP = load(fullfile(param.bp_path,param.bp_file));
if strcmp(param.nbwb,'wb')
    BP.bp_y = interp1(BP.freq_bp,BP.bp,freq_y);
end


%% Fish scattering model
FISH_MODEL = load(fullfile(param.fish_path,param.fish_file));
if strcmp(param.nbwb,'wb')
    FISH_MODEL.fbs_y = interp1(FISH_MODEL.freq_fish,FISH_MODEL.fbs_len_angle,freq_y);
    FISH_MODEL.fbs_y(isnan(FISH_MODEL.fbs_y)) = 0;
end


%% Adjust angles from [deg] to [rad]
if ~strcmp(param.scatterer.type,'rayl') && ~strcmp(param.scatterer.type,'point')
    param.scatterer.angle_mean = param.scatterer.angle_mean/180*pi;
    param.scatterer.angle_std = param.scatterer.angle_std/180*pi;
end


%% Frame length parameters
gateL = round(param.gate_len/param.c/dt);
frameL = gateL+2*yL;
% t_frame = (0:frameL-1)*dt;
mid_frame_pt = round((frameL+1)/2);

% resp_temp = zeros(frameL,num_sample);
s = zeros(1,num_sample);


%% Simulation loop
%H_fish = zeros(length(freq_y),sum(N));
tic
parfor iS=1:num_sample  % realization loop

    % Fish scattering response
    switch param.scatterer.type
        case 'point'   % Point scatterer
            H_fish = [];
            for iN=1:length(N)
                H_fish = [H_fish, mix_r(iN)/sqrt(2)*ones(length(freq_y),sum(N))];
            end
            if strcmp(param.nbwb,'wb')
                H_fish = repmat(H_fish,length(freq_y),1);
            end
            
        case 'rayl'    % Rayleigh scatterer
            H_fish = [];
            for iN=1:length(N)
                H_fish = [H_fish;raylrnd(mix_r(iN)/sqrt(2),1,N(iN))];
            end
            if strcmp(param.nbwb,'wb')
                H_fish = repmat(H_fish,length(freq_y),1);
            end
            
        case 'fish'    % Fish-like scatterer
            % randomly generate fish length
            [len,~] = discrete_rnd(param.scatterer.len_bin,param.scatterer.len_dist,sum(N));          % discrete randomd number
            [~,len_idx] = min(abs(repmat(len,1,length(FISH_MODEL.len))-...  % select fish freq response from preloaded model
                                  repmat(FISH_MODEL.len,sum(N),1)),[],2);
                              
            % randomly generate the angle of orientation
            angle = normrnd_truncated(param.scatterer.angle_mean,param.scatterer.angle_std,2,sum(N),[])'; % truncated normal distribution
            %angle = unifrnd(angle_mean-angle_std,angle_mean+angle_std,sum(N),1); % uniform
            [~,angle_idx] = min(abs(repmat(angle,1,length(FISH_MODEL.angle))-...
                                    repmat(FISH_MODEL.angle,sum(N),1)),[],2);
            
            % select bp response from preloaded model
            idx = (len_idx-1)*length(FISH_MODEL.angle)+angle_idx;
            if strcmp(param.nbwb,'wb')       % broadband fish response
                H_fish = FISH_MODEL.fbs_y(:,idx);
                H_fish(1,:) = 0;  % set freq=0 component=0
            elseif strcmp(param.nbwb,'nb')   % narrowband fish response
                H_fish = interp1(FISH_MODEL.freq_fish,FISH_MODEL.fbs_len_angle,param.nb_freq);
                H_fish = H_fish(idx);
            end

%         case 'prosph'  % Prolate spheroid    % =======NOT FINISHED========
%             % Prolate spheroid high-freq asymptotic solution
%             phase = unifrnd(0,2*pi,N,1);
%             
%             cc = 1;
%             e_ac = 1/10;
%             b1 = cc*e_ac; % length of semi-minor axis
%             if strcmp(param.scatterer.sph_rot_opt,'2D')    % Constrain spheroid rotation in MRA plane
%                 theta_sph = unifrnd(0,2*pi,N,1);  % before 2017/04,
%             elseif strcmp(param.scatterer.sph_rot_opt,'3D')  % theta_sph follow sin(theta_sph) in 3D spherical coord
%                 u = unifrnd(0,1,N,1);
%                 theta_sph = pi/2-acos(u);  % theta_sph calculated from normal incidence
%             end
%             fss = cc/2.*sin(atan(b1./(cc.*tan(theta_sph)))).^2./cos(theta_sph).^2;
%             
%             roughness = raylrnd(ones(N,1)*1/sqrt(2));
%             amp = fss.*roughness;            
%             s = amp.*exp(1i*phase);
% 
%             H_fish = repmat(s,1,length(freq_y))';
%             
%             H_fish = [];
%             for iN=1:length(N)
%                 %amp_tmp = prosph_3D_simulation(ar,N(iN),1,mix_r(iN));  % all angles of orientation
%                 amp_tmp = prosph_amp(ar,N(iN),1,mix_r(iN),...
%                     [angle_mean,angle_std],[len_bin',len_dist']);
%                 H_fish = [H_fish, repmat(amp_tmp,length(freq_y),1)];
%             end
%             %for iN=1:length(N)
%             %    amp_tmp = prosph_3D_simulation(ar,N(iN),1,mix_r(iN));
%             %    H_fish(:,cumsum([0 N(1:iN-1)])+(1:N(iN))) = repmat(amp_tmp,length(freq_y),1);
%             %end

    end

    % Beampattern response
    u = unifrnd(0,1,sum(N),1); 
    theta = acos(u);  % angle in the beam
%     theta = rand_piston_angle(sum(N),bpa)';  % angle in the beam
    [~,ind] = min(abs(repmat(theta,1,length(BP.theta))-...
                      repmat(BP.theta,sum(N),1)),[],2); % pick the right bp idx
    if strcmp(param.nbwb,'wb')
        H_bp = BP.bp_y(:,ind);
    elseif strcmp(param.nbwb,'nb')
%     fac = (2*pi*param.nb_freq/param.c)*0.211*sin(theta);
%     H_bp = (2*besselj(1,fac)./(fac)).^2;
        H_bp = interp1(BP.freq_bp,BP.bp,param.nb_freq);
        H_bp = H_bp(ind);
    end

    % wb: assemble and ifft
    if strcmp(param.nbwb,'wb')
        H_scat = repmat(Rss,1,sum(N)).*H_fish.*H_bp;
        h_scat = ifftshift(ifft([H_scat;flipud(conj(H_scat(2:end,:)))]),1);
    else
        h_scat = y.'*(H_bp.*H_fish);
    end

    % Delay (location of fish)
    % need to do in the time domain since phase variation > 2*pi
    time = round(rand(sum(N),1)*gateL+yL);
    
    % Time domain impulse summation
    resp_temp = zeros(frameL,sum(N));
    for iN=1:sum(N)
        resp_temp(:,iN) = [zeros(time(iN)-yHalfL,1);h_scat(:,iN);...
                           zeros(frameL-time(iN)-yHalfL+1,1)];
    end
    resp_temp = sum(resp_temp,2);
    if ~isreal(resp_temp)
        disp('Error: Invalid time series with imaginary part');
    end
    resp_temp = abs(hilbert(resp_temp));

    s(iS) = resp_temp(mid_frame_pt);
    
end % realization loop

disp('time to generate all samples')
toc
