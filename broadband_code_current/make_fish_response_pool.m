% 2013 07 29  Make pool for fish scattering response
%             include both length variation and angle of orientation
%             variation. Only need high frequency part since the
%             frequency range of AirMarLow is 30-70 kHz

addpath /mnt/storage/fish_scat_model
folder = '/mnt/storage/broadband_code_current/fish_info';

%% Set model parameters
freq_fish = logspace(log10(1e3),log10(500e3),100);

shape = 1;     % 1-prolate spheroid
               % 2-swimbladder shape
               % 3-staight cylinder
opt = 1;       % 1-set V00
               % 2-use default values
swb_back = 1;  % 0-curved back
               % 1-straight back
V00 = 13;      % [cm^3]
D = 170;       % [m]

% g&h from Chu
gamma0 = 1.4;
g_swb = 0.0012*(1+0.1*D)^gamma0;              % rho2/rho1
h_swb = 0.22*(1+0.1*D)^((1-gamma0)/2);        % c2/c1

% angle and length bins
fishL = load(['/mnt/storage/broadband_code_current/fish_info/' ...
              'fish_len_dist.mat']);
len = fishL.L_bin;       % [m]
dangle = 0.1;            % [deg]
angle = (-90:dangle:90)/180*pi; % [rad]

fname = sprintf('fish_scat_response_angle%dto%ddeg_dangel%2.2f_len%dto%dcm.mat',...
                angle(1)/pi*180,angle(end)/pi*180,dangle,...
                round(len(1)*100), round(len(end)*100));

%% Make pool
for iL = 1:length(len)
    disp(['fish length = ',num2str(len(iL)*100),' cm']);

    % prolate spheroid shape position vector
    [xg,ar] = swb_depth_compress(len(iL),D,shape,opt,V00);
    if swb_back==0
        disp('curved back');
        rpos = [xg',ones(size(xg'))];
    else
        disp('straight back');
        rpos = [xg',ar'];
    end

    % angle of orientation variation
    parfor iA = 1:length(angle)
        disp(['angle = ',num2str(angle(iA)/pi*180),' deg']);
        warning off all
        fbs_angle(:,iA) = fbs_deformedCyl_general_new(freq_fish,angle(iA),g_swb,h_swb,ar,rpos);
    end
    fbs_len_angle(:,(iL-1)*length(angle)+(1:length(angle))) = fbs_angle;
end

freq_fish = freq_fish;
param.shape = shape;
param.opt = opt;
param.swb_back = swb_back;
param.V00 = V00;
param.D = D;
param.gamma0 = gamma0;
param.g_swb = g_swb;
param.h.swb = h_swb;

save([folder,'/',fname],'freq_fish','fbs_len_angle','len','angle','param');
