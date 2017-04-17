% Compare broadband echo pdf model predictions using Rayleigh scatterer
% and rough prolate spheroid scatterer model for individuals

SAVE_DIR = '/mnt/storage/bb_echopdf_figs/rayl_prosph_echopdf_cmp';

N_all=[10:10:100,150:50:1000];

rayl_DIR = '/mnt/storage/ECHO_STAT/20121116_pdfmodel_w_respenv_5e4';
rayl_fpost = 'gateLen0.5_freqDepBP.mat';
[pm_rayl,pmx_rayl,N_rayl,s_all_rayl] = get_model_curves_kde(rayl_DIR,N_all,1e2,rayl_fpost);

ps_DIR = '/mnt/storage/ECHO_STAT/20121018_bb_smpl5e4_prosph';
ps_fpost = 'gateLdist0.5_freqDepBP.mat';
[pm_ps,pmx_ps,N_ps,s_all_ps] = get_model_curves_kde(ps_DIR,N_all,1e2,ps_fpost);

raylx = logspace(-4,2,100);
rayl = raylpdf(raylx,1/sqrt(2));

% plot
for N=[10,50,100,300,500,1000]
    [~,idx] = min(abs(N_all-N));
    figure;
    loglog(raylx,rayl,'color',[1 1 1]*180/255,'linewidth',1);
    hold on
    loglog(pmx_rayl,pm_rayl(:,idx),'k','linewidth',0.5);
    loglog(pmx_ps,pm_ps(:,idx),'k--','linewidth',1.5);
    axis([1e-3 1e2 1e-5 1e2]);
    legend('Rayleigh','Rayleigh scatterer','Prolate spheroid scatterer');
    title(['N=',num2str(N)],'fontsize',12);
    set(gca,'fontsize',12,'ticklength',[0.02 0.02]);
    xlabel('Normalized echo amplitude','fontsize',12);
    ylabel('PDF','fontsize',12);
    saveas(gcf,[SAVE_DIR,'/rayl_prosph_cmp_N',num2str(N),'.fig'],'fig');
    saveas(gcf,[SAVE_DIR,'/rayl_prosph_cmp_N',num2str(N),'.png'],'png');
end
