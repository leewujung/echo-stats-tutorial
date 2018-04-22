% Code to generate Figure 5 of the echo statistics tutorial
% This figure plots the Rayleigh PDF and the associated CDF and PFA in both
% lin-lin and log-log scales 
%
% @author: Wu-Jung Lee | leewujung@gmail.com | APL-UW


clear

addpath './util_fcn'
save_base_path = './figs';

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
saveas(gcf,fullfile(save_path,[script_name,'_loglog.fig']),'fig');
saveSameSize(gcf,'file',fullfile(save_path,[script_name,'_loglog.png']),...
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
ll_pos = get(ll,'position');
set(ll,'position',[ll_pos(1),0.55,ll_pos(3),ll_pos(4)])
% title('Rayleigh, lin-lin');
saveas(gcf,fullfile(save_path,[script_name,'_linlin.fig']),'fig');
saveSameSize(gcf,'file',fullfile(save_path,[script_name,'_linlin.png']),...
    'format','png','renderer','painters');
