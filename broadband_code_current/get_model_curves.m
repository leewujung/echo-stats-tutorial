% 2012 05 15  Compile all modeling results into one matrix to facilitate
%             curve-fitting

%{
function [pmx,pm] = get_model_curves(N,sdir,ns,glen,fpre,npt)

tail = sprintf('_sampleN%d_gateLden%2.1f_freqDepBP.mat',...
               ns,glen);
ftail = sprintf('_sampleN%d_gateLdist%2.1f_freqDepBP.mat',...
                ns,glen);

pmx = logspace(-5,100,1000);
pm = zeros(length(pmx),length(N));
for iN=1:length(N)
    fname = [fpre,'N_',num2str(N(iN)),'_r_1',ftail];
    %fname = [fpre,'N',num2str(N(iN)),ftail];
    A = load([sdir,'/',fname]);
    [pm_tmp,pmx_tmp] = findEchoDist_kde(A.s/sqrt(mean(A.s.^2)),npt);
    pm(:,iN) = interp1(pmx_tmp,pm_tmp,pmx);
end
%}

for iN=1:length(N)
    figure;
    loglog(raylx,rayl,'k')
    hold on
    ylim([1e-7 1e2]);
    xlim([1e-3 1e2])
    loglog(pmx_fish,pm_fish(:,iN));
    loglog(pmx_rayl,pm_rayl(:,iN),'r');
    legend('Rayleigh distr','fish','Rayleigh scatterer')
    xlabel('Normalized echo amplitude')
    ylabel('PDF')
    title(['fish and Rayleigh scatterer cmp, N=',num2str(N(iN))]);
    saveas(gcf,['/mnt/storage/ECHO_STAT/20130807_bbechopdf/fish_rayl_cmp_N',num2str(N(iN)),'.fig'],'fig');
    saveas(gcf,['/mnt/storage/ECHO_STAT/20130807_bbechopdf/fish_rayl_cmp_N',num2str(N(iN)),'.png'],'png');
end
