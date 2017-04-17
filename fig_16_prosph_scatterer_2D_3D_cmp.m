% 2017 01 01  Prolate spheroid scatterer with and without beampattern, PDF & PFA
% 2017 04 17  1) Change prolate spheroid orientaiton to 3D instead of 2D
%                2D-rotating within plane containing MRA
%                3D-follow distribution in a spherical coordinate
%             2) Update polar angle generation to use inverse tranform sampling
%             3) Update figure legend, axis labels, and curve style

clear

addpath '~/Dropbox/0_CODE/MATLAB/saveSameSize/'
addpath '~/Dropbox/0_CODE/prolatespheroid/'

% base_path = '~/Desktop/echo_stat_figs';
base_path = '/Volumes/wjlee_apl_2/echo_stat_tutorial/echo_stat_figs/';

% Make save path
str = strsplit(mfilename('fullpath'),'/');
str = str{end};
save_path = fullfile(base_path,str);
if ~exist(save_path,'dir')
    mkdir(save_path);
end
data_path = 'fig_16_prosph_scatterer';

% Set param
X = load('fig_12_pb_ka_ka_num.mat');
ka = X.ka_3deg;
% ka = 2*pi;
sph_rot_opt = '2D';

pingnum_str = '1e7';
pingnum = eval(pingnum_str);
npt = 200;  % number of points for pe kde estimation
N_all = [1,10,100,1000];


% Plot: PDF NO BEAMPATTERN
for iN=1:length(N_all)
    fig = figure;
    xr = logspace(-3,log10(2000),500);  % standard
    rayl = raylpdf(xr,1/sqrt(2));
    loglog(xr,rayl,'k','linewidth',1);
    hold on
    
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_2D_bp0.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(base_path,data_path,simu_file));
    [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
    loglog(x,p_x,'-','linewidth',1);
    
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_3D_bp0.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(base_path,data_path,simu_file));
    [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
    loglog(x,p_x,'-','linewidth',1);
    
    title(sprintf('ka=%2.4f, smplN=%s, N=%d, no bp',...
        ka,pingnum_str,N_all(iN)),...
        'fontsize',18);
    ll = legend('Rayleigh','2D (half-plane)','3D (half-space)',...
        'location','southwest');
    set(ll,'fontsize',18);
    set(gca,'fontsize',16)
    xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
    ylabel('$p_e(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
    xlim([1e-3 1e2]);
    ylim([1e-6 1e3]);
    
    save_fname = sprintf('%s_ka%2.4f_smpl%s_N%04d_pdf_bp0',...
        str,ka,pingnum_str,N_all(iN));
    saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
    saveSameSize_100(fig,'file',[fullfile(save_path,save_fname),'.png'],...
        'format','png');
end



% Plot: PDF WITH BEAMPATTERN
for iN=1:length(N_all)
    fig = figure;
    xr = logspace(-3,log10(2000),500);  % standard
    rayl = raylpdf(xr,1/sqrt(2));
    loglog(xr,rayl,'k','linewidth',2);
    hold on
    
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_2D_bp1.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(base_path,data_path,simu_file));
    [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
    loglog(x,p_x,'-','linewidth',1);
    
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_3D_bp1.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(base_path,data_path,simu_file));
    [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
    loglog(x,p_x,'-','linewidth',1);
    
    title(sprintf('ka=%2.4f, smplN=%s, N=%d, with bp',...
        ka,pingnum_str,N_all(iN)),...
        'fontsize',18);
    ll = legend('Rayleigh','2D (half-plane)','3D (half-space)',...
        'location','southwest');
    set(ll,'fontsize',18);
    set(gca,'fontsize',16)
    xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
    ylabel('$p_e(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
    xlim([1e-3 1e3]);
    ylim([1e-7 1e3]);
    
    save_fname = sprintf('%s_ka%2.4f_smpl%s_N%04d_pdf_bp1',...
        str,ka,pingnum_str,N_all(iN));
    saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
    saveSameSize_100(fig,'file',[fullfile(save_path,save_fname),'.png'],...
        'format','png');
end


% Plot: PFA WITH BEAMPATTERN
for iN=1:length(N_all)
    fig = figure;
    xr = logspace(-3,10,5000);
    rayl = raylpdf(xr,1/sqrt(2));
    cdf_rayl = cumtrapz(xr,rayl);
    pfa_rayl = 1-cdf_rayl;
    loglog(xr,pfa_rayl,'k','linewidth',2);
    hold on
    
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_2D_bp1.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(base_path,data_path,simu_file));
    [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
    cdf_x = cumtrapz(x,p_x);
    pfa_x = 1-cdf_x;
    loglog(x,pfa_x,'-','linewidth',1);
    
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_3D_bp1.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(base_path,data_path,simu_file));
    [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
    cdf_x = cumtrapz(x,p_x);
    pfa_x = 1-cdf_x;
    loglog(x,pfa_x,'-','linewidth',1);
    
    title(sprintf('ka=%2.4f, smplN=%s, N=%d, with bp',...
        ka,pingnum_str,N_all(iN)),...
        'fontsize',18);
    ll = legend('Rayleigh','2D (half-plane)','3D (half-space)',...
        'location','southwest');
    set(ll,'fontsize',18);
    set(gca,'fontsize',16)
    xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
    ylabel('$PFA(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
    xlim([1e-3 1e3]);
    ylim([1e-6 1e3]);
    
    save_fname = sprintf('%s_ka%2.4f_smpl%s_N%04d_pfa_bp1',...
        str,ka,pingnum_str,N_all(iN));
    saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
    saveSameSize_100(fig,'file',[fullfile(save_path,save_fname),'.png'],...
        'format','png');
    
end




% Plot: *RAW* PDF WITH BEAMPATTERN
for iN=1:length(N_all)
    fig = figure;
    xr = logspace(-3,log10(2000),500);  % standard
    rayl = raylpdf(xr,1/sqrt(2));
    loglog(xr,rayl,'k','linewidth',2);
    hold on
    
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_2D_bp1.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(base_path,data_path,simu_file));
    [p_x,x] = findEchoDist_kde(E.env,npt);
    loglog(x,p_x,'-','linewidth',1);
    
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_3D_bp1.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(base_path,data_path,simu_file));
    [p_x,x] = findEchoDist_kde(E.env,npt);
    loglog(x,p_x,'-','linewidth',1);
    
    title(sprintf('ka=%2.4f, smplN=%s, N=%d, with bp',...
        ka,pingnum_str,N_all(iN)),...
        'fontsize',18);
    ll = legend('Rayleigh','2D (half-plane)','3D (half-space)',...
        'location','southwest');
    set(ll,'fontsize',18);
    set(gca,'fontsize',16)
    xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
    ylabel('$p_e(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
    xlim([1e-3 1e1]);
    ylim([1e-7 1e3]);
    
    save_fname = sprintf('%s_ka%2.4f_smpl%s_N%04d_pdf_bp1_raw',...
        str,ka,pingnum_str,N_all(iN));
    saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
    saveSameSize_100(fig,'file',[fullfile(save_path,save_fname),'.png'],...
        'format','png');
end
