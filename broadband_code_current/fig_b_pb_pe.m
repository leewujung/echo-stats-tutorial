% 2012 11 07  Plot b, pb, and pe

SAVE_DIR = '/mnt/storage/bb_echopdf_figs/b_pb_pe';

a = 0.054;
c = 1500;
ka30 = 2*pi*30e3/c*a;
ka50 = 2*pi*50e3/c*a;
ka70 = 2*pi*70e3/c*a;

dtheta = 1e-4*pi;
theta = 0:dtheta:pi/2;
theta(1) = 1e-16;
btheta_30k = (2*besselj(1,ka30*sin(theta))./(ka30*sin(theta))).^2;
btheta_50k = (2*besselj(1,ka50*sin(theta))./(ka50*sin(theta))).^2;
btheta_70k = (2*besselj(1,ka70*sin(theta))./(ka70*sin(theta))).^2;

[b_30k,pb_30k] = calc_pb(ka30,1e-5);
[b_50k,pb_50k] = calc_pb(ka50,1e-5);
[b_70k,pb_70k] = calc_pb(ka70,1e-5);

[e_30k,pe_30k] = calc_pe(1e-5,100,10000,b_30k,pb_30k);
[e_50k,pe_50k] = calc_pe(1e-5,100,10000,b_50k,pb_50k);
[e_70k,pe_70k] = calc_pe(1e-5,100,10000,b_70k,pb_70k);

save('/mnt/storage/bb_echopdf_figs/b_pb_pe/b_pb_pe.mat');

addpath /mnt/storage/modeling_code_current
[b_30k,pb_30k] = normpdfx_xstr(b_30k,pb_30k);
[b_50k,pb_50k] = normpdfx_xstr(b_50k,pb_50k);
[b_70k,pb_70k] = normpdfx_xstr(b_70k,pb_70k);
[e_30k,pe_30k] = normpdfx_xstr(e_30k,pe_30k);
[e_50k,pe_50k] = normpdfx_xstr(e_50k,pe_50k);
[e_70k,pe_70k] = normpdfx_xstr(e_70k,pe_70k);

xrayl = b_30k;
rayl = raylpdf(xrayl,1/sqrt(2));

% b(theta)
figure;
plot(theta/pi*180,20*log10(btheta_30k),'k');
xlabel('{\theta}','fontsize',14);
ylabel('B(\theta) (dB)','fontsize',14);
set(gca,'xtick',0:30:90,'ytick',-80:40:0,'fontsize',12);
axis([0 90 -80 0]);
title('30 kHz','fontsize',14)
saveas(gcf,[SAVE_DIR,'/bthetat_30.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/bthetat_30.png'],'png');

figure;
plot(theta/pi*180,20*log10(btheta_50k),'k');
xlabel('{\theta}','fontsize',14);
ylabel('B(\theta) (dB)','fontsize',14);
set(gca,'xtick',0:30:90,'ytick',-100:50:0,'fontsize',12);
axis([0 90 -100 0]);
title('50 kHz','fontsize',14)
saveas(gcf,[SAVE_DIR,'/bthetat_50.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/bthetat_50.png'],'png');

figure;
plot(theta/pi*180,20*log10(btheta_70k),'k');
xlabel('{\theta}','fontsize',14);
ylabel('B(\theta) (dB)','fontsize',14);
set(gca,'xtick',0:30:90,'ytick',-100:50:0,'fontsize',12);
axis([0 90 -100 0]);
title('70 kHz','fontsize',14)
saveas(gcf,[SAVE_DIR,'/bthetat_70.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/bthetat_70.png'],'png');

% pB(b) & pE(e)
figure;
loglog(xrayl,rayl,'color',[1 1 1]*180/255,'linewidth',2);
hold on
loglog(b_30k,pb_30k,'k','linewidth',1);
loglog(e_30k,pe_30k,'k','linewidth',2);
axis([1e-4 1e2 1e-4 3e3]);
xlabel('Normalized echo amplitude','fontsize',14);
ylabel('PDF','fontsize',14);
set(gca,'fontsize',12);
legend('Rayleigh','P_B(b)','p_E(e)');
title('30 kHz','fontsize',14)
saveas(gcf,[SAVE_DIR,'/pbpe_30.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/pbpe_30.png'],'png');

figure;
loglog(xrayl,rayl,'color',[1 1 1]*180/255,'linewidth',2);
hold on
loglog(b_50k,pb_50k,'k','linewidth',1);
loglog(e_50k,pe_50k,'k','linewidth',2);
axis([1e-4 1e2 1e-4 3e3]);
xlabel('Normalized echo amplitude','fontsize',14);
ylabel('PDF','fontsize',14);
set(gca,'fontsize',12);
legend('Rayleigh','P_B(b)','p_E(e)');
title('50 kHz','fontsize',14)
saveas(gcf,[SAVE_DIR,'/pbpe_50.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/pbpe_50.png'],'png');

figure;
loglog(xrayl,rayl,'color',[1 1 1]*180/255,'linewidth',2);
hold on
loglog(b_70k,pb_70k,'k','linewidth',1);
loglog(e_70k,pe_70k,'k','linewidth',2);
axis([1e-4 1e2 1e-4 3e3]);
xlabel('Normalized echo amplitude','fontsize',14);
ylabel('PDF','fontsize',14);
set(gca,'fontsize',12);
legend('Rayleigh','P_B(b)','p_E(e)');
title('70 kHz','fontsize',14)
saveas(gcf,[SAVE_DIR,'/pbpe_70.fig'],'fig');
saveas(gcf,[SAVE_DIR,'/pbpe_70.png'],'png');



