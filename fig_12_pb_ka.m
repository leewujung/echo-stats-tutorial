% Fig. 12 of echo stat tutorial
% 2016 10 18  first plot
% 2017 04 10  Update curve style and figure legend, axis labels

addpath '~/Dropbox/0_CODE/MATLAB/saveSameSize/'

% save_base_path = '~/Desktop/echo_stat_figs';
save_base_path = '/Volumes/wjlee_apl_2/echo_stat_tutorial/echo_stat_figs/';

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(save_base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

calc_ka_opt = 0;

if calc_ka_opt
    % Find rough ka region
    ka = 6.65:1e-2:6.7;
    theta_deg = 9.5:1e-3:10.5;
    theta = theta_deg/180*pi;
    btheta = (2*besselj(1,ka'*sin(theta))./(ka'*sin(theta))).^2;
    btheta_log = 20*log10(btheta);
    
    figure
    plot(theta/pi*180,btheta_log)
    xlabel('Polar angle (^o)');
    ylabel('2-way beampattern');
    grid
    
    % Full beamwidth 1 deg
    ka = 132.7:1e-5:132.8;
    theta_deg_want = 0.5;
    theta_deg = theta_deg_want-0.005:1e-6:theta_deg_want+0.005;
    theta = theta_deg/180*pi;
    btheta = (2*besselj(1,ka'*sin(theta))./(ka'*sin(theta))).^2;
    btheta_log = 20*log10(btheta);
    
    [~,idx_3db] = min(abs(btheta_log-(-3)),[],2);
    [~,idx_ka] = min(abs(theta_deg(idx_3db)-theta_deg_want));
    
    ka_want = ka(idx_ka);
    btheta_want = (2*besselj(1,ka_want'*sin(theta))./(ka_want'*sin(theta))).^2;
    
    ka_1deg = ka_want;
    
    
    % Full beamwidth 5 deg
    ka = 26.5:1e-5:26.6;
    theta_deg_want = 2.5;
    theta_deg = theta_deg_want-0.005:1e-6:theta_deg_want+0.005;
    theta = theta_deg/180*pi;
    btheta = (2*besselj(1,ka'*sin(theta))./(ka'*sin(theta))).^2;
    btheta_log = 20*log10(btheta);
    
    [~,idx_3db] = min(abs(btheta_log-(-3)),[],2);
    [~,idx_ka] = min(abs(theta_deg(idx_3db)-theta_deg_want));
    
    ka_want = ka(idx_ka);
    btheta_want = (2*besselj(1,ka_want'*sin(theta))./(ka_want'*sin(theta))).^2;
    
    ka_5deg = ka_want;
    
    
    % Full beamwidth 10 deg
    ka = 13.2:1e-5:13.3;
    theta_deg_want = 5;
    theta_deg = theta_deg_want-0.005:1e-6:theta_deg_want+0.005;
    theta = theta_deg/180*pi;
    btheta = (2*besselj(1,ka'*sin(theta))./(ka'*sin(theta))).^2;
    btheta_log = 20*log10(btheta);
    
    [~,idx_3db] = min(abs(btheta_log-(-3)),[],2);
    [~,idx_ka] = min(abs(theta_deg(idx_3db)-theta_deg_want));
    
    ka_want = ka(idx_ka);
    btheta_want = (2*besselj(1,ka_want'*sin(theta))./(ka_want'*sin(theta))).^2;
    
    ka_10deg = ka_want;
    
    
    % Full beamwidth 3 deg
    ka = 44.2:1e-5:44.3;
    theta_deg_want = 1.5;
    theta_deg = theta_deg_want-0.005:1e-6:theta_deg_want+0.005;
    theta = theta_deg/180*pi;
    btheta = (2*besselj(1,ka'*sin(theta))./(ka'*sin(theta))).^2;
    btheta_log = 20*log10(btheta);
    
    [~,idx_3db] = min(abs(btheta_log-(-3)),[],2);
    [~,idx_ka] = min(abs(theta_deg(idx_3db)-theta_deg_want));
    
    ka_want = ka(idx_ka);
    btheta_want = (2*besselj(1,ka_want'*sin(theta))./(ka_want'*sin(theta))).^2;
    
    ka_3deg = ka_want;
    
    
     % Full beamwidth 20 deg
    ka = 6.65:1e-5:6.7;
    theta_deg_want = 10;
    theta_deg = theta_deg_want-0.005:1e-6:theta_deg_want+0.005;
    theta = theta_deg/180*pi;
    btheta = (2*besselj(1,ka'*sin(theta))./(ka'*sin(theta))).^2;
    btheta_log = 20*log10(btheta);
    
    [~,idx_3db] = min(abs(btheta_log-(-3)),[],2);
    [~,idx_ka] = min(abs(theta_deg(idx_3db)-theta_deg_want));
    
    ka_want = ka(idx_ka);
    btheta_want = (2*besselj(1,ka_want'*sin(theta))./(ka_want'*sin(theta))).^2;
    
    ka_20deg = ka_want;
    
    
    % Save ka numbers
    save(fullfile(save_path,[script_name,'_ka_num.mat']),...
        'ka_*deg');
    
end


% Plot pb
if 0
    b_num = 5e4;
    b_start_log = -7;
    b_end_log = 0;
    [b_1deg,pb_1deg] = calc_pb_log(ka_1deg,b_start_log,b_end_log,b_num);
    [b_3deg,pb_3deg] = calc_pb_log(ka_3deg,b_start_log,b_end_log,b_num);
    [b_5deg,pb_5deg] = calc_pb_log(ka_5deg,b_start_log,b_end_log,b_num);
    [b_10deg,pb_10deg] = calc_pb_log(ka_10deg,b_start_log,b_end_log,b_num);
    % [b_20deg,pb_20deg] = calc_pb_log(ka_20deg,b_start_log,b_end_log,b_num);
    
    save(fullfile(save_path,[script_name,'_b_pb.mat']),...
        'b_*deg','pb_*deg');
end

load(fullfile(save_path,[script_name,'_b_pb.mat']));

% Figures
figure
loglog(b_10deg,pb_10deg)
hold on
loglog(b_5deg,pb_5deg)
loglog(b_3deg,pb_3deg)
loglog(b_1deg,pb_1deg)
axis([1e-7 1e0 1e-4 1e8])
xlabel('Echo amplitude');
ylabel('PDF');
legend('20^o','10^o','5^o','3^o','1^o','location','southwest');
title('P_b(b), different full beamwidth');
saveas(gcf,fullfile(save_path,[script_name,'_b_pb_all.fig']),'fig');
saveSameSize(gcf,'file',fullfile(save_path,[script_name,'_b_pb_all.png']),...
    'format','png','renderer','painters');

figure
loglog(b_10deg,pb_10deg,'k:','linewidth',1.5)
hold on
loglog(b_5deg,pb_5deg,'k--')
loglog(b_3deg,pb_3deg,'k-.')
loglog(b_1deg,pb_1deg,'k');
axis([1e-7 1e0 1e-4 1e8])
xlabel('Echo amplitude');
ylabel('PDF');
legend('20^o','10^o','5^o','3^o','1^o','location','southwest');
title('P_b(b), different full beamwidth');
saveas(gcf,fullfile(save_path,[script_name,'_b_pb_all_bw.fig']),'fig');
saveSameSize(gcf,'file',fullfile(save_path,[script_name,'_b_pb_all_bw.png']),...
    'format','png','renderer','painters');

% Individual Pb plots
figure
loglog(b_1deg,pb_1deg,'k','linewidth',1)
axis([1e-7 2e0 1e-4 1e8])
set(gca,'fontsize',16)
ylabel('$p_b(b)$','Interpreter','LaTex','fontsize',24);
xlabel('$b$','Interpreter','LaTex','fontsize',24);
% title('P_b(b), full beamwidth = 1^o');
saveas(gcf,fullfile(save_path,[script_name,'_b_pb_1deg.fig']),'fig');
saveSameSize(gcf,'file',fullfile(save_path,[script_name,'_b_pb_1deg.png']),...
    'format','png','renderer','painters');

figure
loglog(b_3deg,pb_3deg,'k','linewidth',1)
axis([1e-7 2e0 1e-4 1e8])
set(gca,'fontsize',16)
ylabel('$p_b(b)$','Interpreter','LaTex','fontsize',24);
xlabel('$b$','Interpreter','LaTex','fontsize',24);
% title('P_b(b), full beamwidth = 3^o');
saveas(gcf,fullfile(save_path,[script_name,'_b_pb_3deg.fig']),'fig');
saveSameSize(gcf,'file',fullfile(save_path,[script_name,'_b_pb_3deg.png']),...
    'format','png','renderer','painters');

figure
loglog(b_5deg,pb_5deg,'k','linewidth',1)
axis([1e-7 2e0 1e-4 1e8])
set(gca,'fontsize',16)
ylabel('$p_b(b)$','Interpreter','LaTex','fontsize',24);
xlabel('$b$','Interpreter','LaTex','fontsize',24);
% title('P_b(b), full beamwidth = 5^o');
saveas(gcf,fullfile(save_path,[script_name,'_b_pb_5deg.fig']),'fig');
saveSameSize(gcf,'file',fullfile(save_path,[script_name,'_b_pb_5deg.png']),...
    'format','png','renderer','painters');

figure
loglog(b_10deg,pb_10deg,'k','linewidth',1)
axis([1e-7 2e0 1e-4 1e8])
set(gca,'fontsize',16)
ylabel('$p_b(b)$','Interpreter','LaTex','fontsize',24);
xlabel('$b$','Interpreter','LaTex','fontsize',24);
% title('P_b(b), full beamwidth = 10^o');
saveas(gcf,fullfile(save_path,[script_name,'_b_pb_10deg.fig']),'fig');
saveSameSize(gcf,'file',fullfile(save_path,[script_name,'_b_pb_10deg.png']),...
    'format','png','renderer','painters');



% For verification of ka numbers
theta = (0:1e-3:12)/180*pi;
btheta_1deg = (2*besselj(1,ka_1deg*sin(theta))./(ka_1deg*sin(theta))).^2;
btheta_3deg = (2*besselj(1,ka_3deg*sin(theta))./(ka_3deg*sin(theta))).^2;
btheta_5deg = (2*besselj(1,ka_5deg*sin(theta))./(ka_5deg*sin(theta))).^2;
btheta_10deg = (2*besselj(1,ka_10deg*sin(theta))./(ka_10deg*sin(theta))).^2;
% btheta_20deg = (2*besselj(1,ka_20deg*sin(theta))./(ka_20deg*sin(theta))).^2;

figure
plot(theta/pi*180,20*log10(btheta_20deg),'k','linewidth',1.5)
hold on
plot(theta/pi*180,20*log10(btheta_10deg),'k:','linewidth',1.5)
plot(theta/pi*180,20*log10(btheta_5deg),'k--');
plot(theta/pi*180,20*log10(btheta_3deg),'k-.');
plot(theta/pi*180,20*log10(btheta_1deg),'k');
xlabel('Polar angle (^o)')
ylabel('2-way beampattern (dB)');
title('Beampattern for different ka/beamwidth')
ylim([-5 0])
grid
legend('20^o beam','10^o beam','5^o beam','3^o beam','1^o beam');
saveas(gcf,fullfile(save_path,[script_name,'_bp_all.fig']),'fig');
saveSameSize(gcf,'file',fullfile(save_path,[script_name,'_bp_all.png']),...
    'format','png','renderer','painters');



