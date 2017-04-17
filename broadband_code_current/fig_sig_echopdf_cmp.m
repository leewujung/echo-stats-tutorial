% 2012 11 08  Plot echo pdf models resulted from different signals

MD_sys = '/mnt/storage/ECHO_STAT/20121116_pdfmodel_w_respenv_5e4';
MD_nosys = '/mnt/storage/ECHO_STAT/20121127_tx_nosys_bbmodel';
%MD_nosys = '/mnt/storage/ECHO_STAT/20121128_tx_nosys_bbmodel_2ndset';
MD_hann = '/mnt/storage/ECHO_STAT/20121127_hannwin_hf_lf';
MD_midh = '/mnt/storage/ECHO_STAT/20121129_midhann_width';
SD = '/mnt/storage/bb_echopdf_figs/sig_echopdf_cmp';
fpost = '_gateLen0.5_freqDepBP.mat';
npt = 1e2;

save_opt = 0;
addpath /mnt/storage/analysis_code_current

xrayl = logspace(-5,2,1e2);
rayl = raylpdf(xrayl,1/sqrt(2));

N = [10,50];
for iN=1:length(N)
    idealtx = load([MD_nosys,'/tx_nosys_N',num2str(N(iN)),'_sampleN50000',fpost]);
    actualtx = load([MD_sys,'/tx_sys_N',num2str(N(iN)),'_sampleN50000',fpost]);
    hihann = load([MD_hann,'/hann_hf_N',num2str(N(iN)),'_sampleN50000',fpost]);
    lohann = load([MD_hann,'/hann_lf_N',num2str(N(iN)),'_sampleN50000',fpost]);
    mhannn = load([MD_midh,'/midhann_narrow_N',num2str(N(iN)),'_sampleN50000',fpost]);
    mhannw = load([MD_midh,'/midhann_wide_N',num2str(N(iN)),'_sampleN50000',fpost]);

    s{iN}(1,:) = idealtx.s(1:5e4);
    s{iN}(2,:) = actualtx.s(1:5e4);
    s{iN}(3,:) = hihann.s(1:5e4);
    s{iN}(4,:) = lohann.s(1:5e4);
    s{iN}(5,:) = mhannn.s(1:5e4);
    s{iN}(6,:) = mhannw.s(1:5e4);
    s_rms{iN} = mean(sqrt(s{iN}.^2),2);
    s{iN} = s{iN}./repmat(s_rms{iN},1,5e4);
    
    for ii=1:6
    [dens{iN}(ii,:),x{iN}(ii,:),bw(iN,1)] = findEchoDist_kde(s{iN}(ii,:),npt);
    end
    
    % ideal vs actual tx
    figure;
    hia = loglog(xrayl,rayl,'color',[1 1 1]*180/255,'linewidth',2);
    hold on
    loglog(x{iN}(1,:),dens{iN}(1,:),'k--','linewidth',3); % ideal
    loglog(x{iN}(2,:),dens{iN}(2,:),'k-','linewidth',1);  % actual
    axis([5e-3 1e2 1e-5 1e1]);
    xlabel('Normalized echo amplitude','fontsize',14);
    ylabel('PDF','fontsize',14);
    legend('Rayleigh','Ideal tx','Actual tx');
    set(gca,'fontsize',12,'ticklength',[1 1]*0.02);
    title(['N=',num2str(N(iN))]);

    % hihann vs lohann
    figure;
    hhl = loglog(xrayl,rayl,'color',[1 1 1]*180/255,'linewidth',2);
    hold on
    loglog(x{iN}(4,:),dens{iN}(4,:),'k','linewidth',1);   % low
    loglog(x{iN}(3,:),dens{iN}(3,:),'k--','linewidth',3); % high
    axis([5e-3 1e2 1e-5 1e1]);
    xlabel('Normalized echo amplitude','fontsize',14);
    ylabel('PDF','fontsize',14);
    legend('Rayleigh','Low Hann win','High Hann win');
    set(gca,'fontsize',12,'ticklength',[1 1]*0.02);
    title(['N=',num2str(N(iN))]);

    
    % narrow vs wide hann
    hwn = figure;
    loglog(xrayl,rayl,'color',[1 1 1]*180/255,'linewidth',2);
    hold on
    loglog(x{iN}(5,:),dens{iN}(5,:),'k','linewidth',1);   % narrow
    loglog(x{iN}(6,:),dens{iN}(6,:),'k--','linewidth',3); % wide
    axis([5e-3 1e2 1e-5 1e1]);
    xlabel('Normalized echo amplitude','fontsize',14);
    ylabel('PDF','fontsize',14);
    legend('Rayleigh','Narrow Hann win','Wide Hann win');
    set(gca,'fontsize',12,'ticklength',[1 1]*0.02);
    title(['N=',num2str(N(iN))]);

    if save_opt==1
    saveas(hia,[SD,'/ideal_actaul_N',num2str(N(iN)),'.fig']);
    saveas(hia,[SD,'/ideal_actaul_N',num2str(N(iN)),'.png']);
    saveas(hhl,[SD,'/hihann_lohann_N',num2str(N(iN)),'.fig']);
    saveas(hhl,[SD,'/hihann_lohann_N',num2str(N(iN)),'.png']);
    saveas(hwn,[SD,'/narrow_wide_hann_N',num2str(N(iN)),'.fig']);
    saveas(hwn,[SD,'/narrow_wide_hann_N',num2str(N(iN)),'.png']);
    end
    
end

