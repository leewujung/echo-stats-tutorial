% 2014 06 26  Plot broadband echo pdfs with and without beampattern effects

addpath /mnt/storage/analysis_code_current

N=[10,50,300];
ns = 5e4;
glen = 0.5;
model_dir_rayleigh_withbp = '/mnt/storage/broadband_echopdf_models/rayleigh_scatterers';
model_dir_rayleigh_nobp = '/mnt/storage/broadband_echopdf_models/rayleigh_scatterers_nobp';

pmx = logspace(-5,100,1000);
rayl = raylpdf(pmx,1/sqrt(2));

npt = 200;
pm_withbp = gather_bbechopdf_curves(N,model_dir_rayleigh_withbp,ns,glen,'rayleigh',npt,pmx);
pm_nobp = gather_bbechopdf_curves(N,model_dir_rayleigh_nobp,ns,glen,'rayleigh_nobp',npt,pmx);

save_folder = '/mnt/storage/broadband_echopdf_ms_figs/fig_with_without_bp';
if ~exist(save_folder,'dir')
    mkdir(save_folder);
end

% FIGURE: rayleigh and fish scatterer comparison
for iN=1:length(N)
    figure;
    loglog(pmx,rayl,'color',[0.5,0.5,0.5]);
    hold on
    loglog(pmx,pm_nobp(:,iN),'b');
    loglog(pmx,pm_withbp(:,iN),'r');
    axis([1e-2 5e1 1e-6 5e1])
    title(sprintf('N=%d',N(iN)));
    xlabel('Normalized echo amplitude');
    ylabel('PDF');
    legend('Rayleigh distr','Without bp','With bp');
    saveas(gcf,sprintf('%s/rayleigh_scatterers_with_without_bp_N%d.fig',...
                       save_folder,N(iN)),'fig');
    print(sprintf('%s/rayleigh_scatterers_with_without_bp_N%d.png',...
                  save_folder,N(iN)),'-dpng');
    
end



