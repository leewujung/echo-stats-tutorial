% Fig. 5 of echo stat tutorial
% 2016 10 18  first plot
% 2017 04 10  Update curve style and figure legend, axis labels

clear

addpath '~/Dropbox/0_CODE'/MATLAB/saveSameSize/
save_base_path = '/Volumes/wjlee_apl_2/echo_stat_tutorial/echo_stat_figs/';

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(save_base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% log-log plot
x_log = logspace(-3,2,120);
rayl_pdf_log = raylpdf(x_log,1/sqrt(2));
rayl_cdf_log = raylcdf(x_log,1/sqrt(2));

figure;
loglog(x_log,rayl_pdf_log,'k','linewidth',2.5);
hold on
loglog(x_log,rayl_cdf_log,'k--','linewidth',1);
loglog(x_log,1-rayl_cdf_log,'k-','linewidth',1);
set(gca,'fontsize',16)
axis([1e-2 1e1 1e-6 1e3])
xlabel('$x/<x^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
ylabel('$\mathcal{F}(x/<x^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
% title('Rayleigh, log-log');
saveas(gcf,fullfile(save_path,[script_name,'_linlin.fig']),'fig');
saveSameSize(gcf,'file',fullfile(save_path,[script_name,'_linlin.png']),...
    'format','png','renderer','painters');

% lin-lin plot
x_lin = linspace(0,10,1000);
rayl_pdf_lin = raylpdf(x_lin,1/sqrt(2));
rayl_cdf_lin = raylcdf(x_lin,1/sqrt(2));

figure;
plot(x_lin,rayl_pdf_lin,'k','linewidth',2.5);
hold on
loglog(x_lin,rayl_cdf_lin,'k--','linewidth',1);
loglog(x_lin,1-rayl_cdf_lin,'k-','linewidth',1);
set(gca,'fontsize',16)
axis([0 3 0 1.15])
xlabel('$x/<x^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
ylabel('$\mathcal{F}(x/<x^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
ll = legend('PDF','CDF','PFA');
set(ll,'fontsize',18)
% title('Rayleigh, lin-lin');
saveas(gcf,fullfile(save_path,[script_name,'_loglog.fig']),'fig');
saveSameSize(gcf,'file',fullfile(save_path,[script_name,'_loglog.png']),...
    'format','png','renderer','painters');
