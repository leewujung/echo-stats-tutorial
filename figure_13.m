% Code to generate Figure 13 of the echo statistics tutorial.
% This figure shows beampattern PDF associated with a circular aperture for
% 3D and 2D distributions of scatterers 
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

pingnum_str = '1e8';
pingnum = eval(pingnum_str);

N_all = 1;
v_rayl = 1/sqrt(2);

% Check if figure_15 results are already in place
% If not ask users to run it first
if ~exist('./figs/figure_15/pnum_1e8_ka44.2511_N0001_bp1.mat','file')
    disp(['Run figure_15 first to obtain necessary simulation for ' ...
          'plotting this figure!'])
    return
end


% Set operation
mc_opt = 0;  % 0 - do not re-generate realizations
             % 1 - re-generate all realizations

% Monte Carlo simulation
if mc_opt 
    for iN=1:length(N_all)
        Ns = N_all(iN);

        param.N = Ns;
        param.ka = ka;

        env_bp1_2D = zeros(1,pingnum);
        parfor iP = 1:pingnum
            phase = rand(1,Ns)*2*pi;
            amp = ones(1,Ns);
            s = amp.*exp(1i*phase);

            % position in 2D beam
            theta_2D = rand(1,Ns)*pi/2;
            b_bp1_2D = (2*besselj(1,ka*sin(theta_2D))./ ...
                        (ka*sin(theta_2D))).^2;

            % E=SB
            e_bp1_2D = s.*b_bp1_2D;
            env_bp1_2D(iP) = abs(sum(e_bp1_2D));

        end % pingnum

        env = env_bp1_2D;
        file_save = sprintf('pnum_%s_ka%2.4f_N%04d_bp1_2D.mat',...
            pingnum_str,ka,Ns);
        save([save_path,'/',file_save],'env','param');
    end
end



% Plot: PDF WITH BEAMPATTERN
save_path_3D = 'figure_15';    % use simulation results from (Fig. 15 point scatterers)
% save_path = 'new_fig_13_2D_3D_bp';
fig = figure;
xr = logspace(-3,log10(2000),500);  % standard
rayl = raylpdf(xr,1/sqrt(2));

for iN=1:length(N_all)    
    simu_file_3D = sprintf('pnum_%s_ka%2.4f_N%04d_bp1.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(base_path,save_path_3D,simu_file_3D));
    [x_3D,p_x_3D] = findEchoDist(E.env,1200);
%     [x_3D,p_x_3D] = findEchoDist(E.env/sqrt(mean(E.env.^2)),1200);
%     cdf_x_3D = cumtrapz(x_3D,p_x_3D);
%     pfa_x_3D = 1-cdf_x_3D;
%     loglog(x_3D,pfa_x_3D,'-','linewidth',1);
    loglog(x_3D,p_x_3D,'k-','linewidth',2);   
    clear E

    hold on

    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_bp1_2D.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(save_path,simu_file));
    [x_2D,p_x_2D] = findEchoDist(E.env,1200);
%     [x_2D,p_x_2D] = findEchoDist(E.env/sqrt(mean(E.env.^2)),1200);
%     cdf_x_2D = cumtrapz(x_2D,p_x_2D);
%     pfa_x_2D = 1-cdf_x_2D;
%     loglog(x_2D,pfa_x_2D,'-','linewidth',1);
    loglog(x_2D,p_x_2D,'k--');
    clear E
    
end
set(gca,'fontsize',16)
ylabel('$p_b(b)$','Interpreter','LaTex','fontsize',24);
xlabel('$b$','Interpreter','LaTex','fontsize',24);
% title(sprintf('ka=%2.4f, smplN=%s, with bp',...
%     ka,pingnum_str),...
%     'fontsize',18);
ll = legend('Half space','Half plane');
set(ll,'fontsize',18);
xlim([1e-6 2]);
ylim([1e-4 1e5]);

save_fname = sprintf('%s_smpl%s_ka%2.4f_pdf',...
    str,pingnum_str,ka);
saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
saveSameSize(fig,'file',[fullfile(save_path,save_fname),'.png'],...
    'format','png');



% 
% % Plot: PFA WITH BEAMPATTERN
% fig = figure;
% xr = logspace(-3,log10(2000),500);  % standard
% rayl = raylpdf(xr,1/sqrt(2));
% loglog(xr,rayl,'k','linewidth',1);
% hold on
% for iN=1:length(N_all)    
%     simu_file_3D = sprintf('pnum_%s_ka%2.4f_N%04d_bp1.mat',...
%         pingnum_str,ka,N_all(iN));
%     E = load(fullfile(base_path,save_path_3D,simu_file_3D));
%     [x_3D,p_x_3D] = findEchoDist(E.env/sqrt(mean(E.env.^2)),600);
% %     [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
%     cdf_x_3D = cumtrapz(x_3D,p_x_3D);
%     pfa_x_3D = 1-cdf_x_3D;
%     loglog(x_3D,pfa_x_3D,'-','linewidth',1);
%     clear E
%     simu_file_2D = sprintf('pnum_%s_ka%2.4f_N%04d_bp1.mat',...
%         pingnum_str,ka,N_all(iN));
%     E = load(fullfile(base_path,save_path,simu_file));
%     [x_2D,p_x_2D] = findEchoDist(E.env/sqrt(mean(E.env.^2)),600);
% %     [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
%     cdf_x_2D = cumtrapz(x_2D,p_x_2D);
%     pfa_x_2D = 1-cdf_x_2D;
%     pfa_x_2D(pfa_x_2D==0) = NaN;
%     loglog(x_2D,pfa_x_2D,'-','linewidth',1);
%     clear E
% end
% xlabel('Normalized echo amplitude','fontsize',16);
% ylabel('PFA','fontsize',16);
% title(sprintf('ka=%2.4f, smplN=%s, with bp',...
%     ka,pingnum_str),...
%     'fontsize',18);
% ll = legend('Rayleigh','2D beampattern','3D beampattern');
% set(ll,'fontsize',18);
% set(gca,'fontsize',14)
% xlim([1e-3 1e2]);
% ylim([1e-6 1e3]);
% 
% save_fname = sprintf('%s_smpl%s_ka%2.4f_pfa',...
%     str,pingnum_str,ka);
% saveas(fig,[fullfile(base_path,save_path,save_fname),'.fig'],'fig');
% saveSameSize_100(fig,'file',[fullfile(base_path,save_path,save_fname),'.png'],...
%     'format','png');
