% Code to generate Figure 17 of the echo statistics tutorial.
%
% This code plots the PDF and PFA of the echo magnitude due to N identical
% randomly rough, randomly oriented prolate spheroids that are randomly
% distributed in the sensor beam (in 17a, the beam is omnidirectional).
% 3D distribution of scatterers.
%
% Author: Wu-Jung Lee | leewujung@gmail.com | APL-UW


clear
addpath './util_fcn'
base_path = './figs';

% Make save path
str = strsplit(mfilename('fullpath'),'/');
str = str{end};
save_path = fullfile(base_path,str);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Set param
X = load('./figs/figure_12/figure_12_ka_num.mat');
ka = X.ka_3deg;

sph_rot_opt = '3D';

pingnum_str = '1e7';
pingnum = eval(pingnum_str);

npt = 200;  % number of points for pe kde estimation

v_rayl = 1/sqrt(2);
N_all = [1,2,5,10,100,1000];

% Set operation
mc_opt = 1;  % 0 - do not re-generate realizations
             % 1 - re-generate all realizations

% Monte Carlo simulation
if mc_opt
    for iN=1:length(N_all)
        Ns = N_all(iN);
        fprintf('Ns = %d\n',Ns);

        param.N = Ns;
        param.ka = ka;

        parfor iP = 1:pingnum
            phase = unifrnd(0,2*pi,Ns,1);

            cc = 1;
            e_ac = 1/10;
            b1 = cc*e_ac; % length of semi-minor axis
            if strcmp(sph_rot_opt,'2D')    % Constrain spheroid rotation in MRA plane
                theta_sph = unifrnd(0,2*pi,Ns,1);  % before 2017/04,
            elseif strcmp(sph_rot_opt,'3D')  % theta_sph follow sin(theta_sph) in 3D spherical coord
                u = unifrnd(0,1,Ns,1);
                theta_sph = pi/2-acos(u);  % theta_sph calculated from normal incidence
            end
            fss = cc/2.*sin(atan(b1./(cc.*tan(theta_sph)))).^2./cos(theta_sph).^2;
            roughness = raylrnd(ones(Ns,1)*1/sqrt(2));
            amp = fss.*roughness;

            s = amp.*exp(1i*phase);

            % position in the beam
            u = unifrnd(0,1,Ns,1);
            theta = acos(u);

            % beampattern modulation
            b_bp1 = (2*besselj(1,ka*sin(theta))./(ka*sin(theta))).^2;
            b_bp0 = 1;

            % E=SB
            e_bp0 = s.*b_bp0;
            e_bp1 = s.*b_bp1;

            env_bp0(iP) = abs(sum(e_bp0));
            env_bp1(iP) = abs(sum(e_bp1));

        end % pingnum loop

        env = env_bp0;
        file_save = sprintf('pnum_%s_ka%2.4f_N%04d_%s_bp0.mat',...
            pingnum_str,ka,Ns,sph_rot_opt);
        save([save_path,'/',file_save],'env','param');

        env = env_bp1;
        file_save = sprintf('pnum_%s_ka%2.4f_N%04d_%s_bp1.mat',...
            pingnum_str,ka,Ns,sph_rot_opt);
        save([save_path,'/',file_save],'env','param');
    end
end



% Plot: PDF NO BEAMPATTERN
N_all = [1,2,5,10];

fig = figure;
xr = logspace(-3,log10(2000),500);  % standard
rayl = raylpdf(xr,1/sqrt(2));
loglog(xr,rayl,'k','linewidth',2);
hold on
for iN=1:length(N_all)
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_%s_bp0.mat',...
        pingnum_str,ka,N_all(iN),sph_rot_opt);
    E = load(fullfile(save_path,simu_file));
    %[x,p_x] = findEchoDist(E.env/sqrt(mean(E.env.^2)),npt);
    [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
    switch iN
        case 1
            loglog(x,p_x,'r-','linewidth',2);
        case 2
            loglog(x,p_x,'g-','linewidth',2);
        case 3
            loglog(x,p_x,'b-','linewidth',2);
        case 4
            loglog(x,p_x,'b-','linewidth',1);
    end
end
% title(sprintf('ka=%2.4f, smplN=%s, no bp',...
%     ka,pingnum_str),...
%     'fontsize',18);
ll = legend('Rayleigh','N=1','N=2','N=5','N=10',...
    'location','southwest');
set(ll,'fontsize',18);
set(gca,'fontsize',16)
xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
ylabel('$p_e(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
xlim([1e-3 1e2]);
ylim([1e-6 1e3]);

save_fname = sprintf('%s_ka%2.4f_smpl%s_%s_pdf_bp0',...
    str,ka,pingnum_str,sph_rot_opt);
saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
saveSameSize(fig,'file',[fullfile(save_path,save_fname),'.png'],...
    'format','png');


% Plot: PDF WITH BEAMPATTERN
N_all = [1,10,100,1000];
fig = figure;
xr = logspace(-3,log10(2000),500);  % standard
rayl = raylpdf(xr,1/sqrt(2));
loglog(xr,rayl,'k','linewidth',2);
hold on
for iN=1:length(N_all)

    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_%s_bp1.mat',...
        pingnum_str,ka,N_all(iN),sph_rot_opt);
    E = load(fullfile(save_path,simu_file));
    %[x,p_x] = findEchoDist(E.env/sqrt(mean(E.env.^2)),npt);
    [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
    switch iN
        case 1
            loglog(x,p_x,'r-','linewidth',2);
        case 2
            loglog(x,p_x,'g-','linewidth',2);
        case 3
            loglog(x,p_x,'b-','linewidth',2);
        case 4
            loglog(x,p_x,'b-','linewidth',1);
    end
end
% title(sprintf('ka=%2.4f, smplN=%s, with bp',...
%     ka,pingnum_str),...
%     'fontsize',18);
ll = legend('Rayleigh','N=1 (0.00375)','N=10 (0.0375)',...
    'N=100 (0.375)','N=1000 (3.75)',...
    'location','southwest');
set(ll,'fontsize',18);
set(gca,'fontsize',16)
xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
ylabel('$p_e(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
xlim([1e-3 1e3]);
ylim([1e-7 1e3]);

save_fname = sprintf('%s_ka%2.4f_smpl%s_%s_pdf_bp1',...
    str,ka,pingnum_str,sph_rot_opt);
saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
saveSameSize(fig,'file',[fullfile(save_path,save_fname),'.png'],...
    'format','png');


% Plot: PFA WITH BEAMPATTERN
N_all = [1,10,100,1000];
fig = figure;
xr = logspace(-3,10,5000);
rayl = raylpdf(xr,1/sqrt(2));
cdf_rayl = cumtrapz(xr,rayl);
pfa_rayl = 1-cdf_rayl;
loglog(xr,pfa_rayl,'k','linewidth',2);
hold on
for iN=1:length(N_all)
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_%s_bp1.mat',...
        pingnum_str,ka,N_all(iN),sph_rot_opt);
    E = load(fullfile(save_path,simu_file));
    %[x,p_x] = findEchoDist(E.env/sqrt(mean(E.env.^2)),npt);
    [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
    cdf_x = cumtrapz(x,p_x);
    pfa_x = 1-cdf_x;
    switch iN
        case 1
            loglog(x,pfa_x,'r-','linewidth',2);
        case 2
            loglog(x,pfa_x,'g-','linewidth',2);
        case 3
            loglog(x,pfa_x,'b-','linewidth',2);
        case 4
            loglog(x,pfa_x,'b-','linewidth',1);
    end
    clear E
end
% title(sprintf('ka=%2.4f, smplN=%s, with bp',...
%     ka,pingnum_str),...
%     'fontsize',18);
ll = legend('Rayleigh','N=1 (0.00375)','N=10 (0.0375)',...
    'N=100 (0.375)','N=1000 (3.75)',...
    'location','southwest');
set(ll,'fontsize',18);
set(gca,'fontsize',16)
xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
ylabel('$PFA(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
xlim([1e-3 1e3]);
ylim([1e-6 1e1]);

save_fname = sprintf('%s_ka%2.4f_smpl%s_%s_pfa_bp1',...
    str,ka,pingnum_str,sph_rot_opt);
saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
saveSameSize(fig,'file',[fullfile(save_path,save_fname),'.png'],...
    'format','png');
