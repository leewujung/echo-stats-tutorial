% 2014 06 25  Compare the broadband echo pdf models calculated
%             using three types of scatterers
% 2014 11 23  Compare the broadband echo pdf models calculated
%             using Rayleigh scatterers and prolate spheroid
%             scatterers, with angle of orientation constraints

% Get bbechopdf model for plotting the pool
%N=[10,20,50,100,300,1000]';
N=[10,20,50,100,300,1000,1500,2000,2500,3000]';
ns = 5e4;
glen = 0.5;
model_dir_prosph = '/mnt/storage/broadband_echopdf_models/prosph_0.05_allAngleDistr_singleLen_all';
model_dir_rayleigh = '/mnt/storage/broadband_echopdf_models/rayleigh_scatterers';
model_dir_fish = '/mnt/storage/broadband_echopdf_models/fish_defaultLenDistr_defaultAngleDistr';

npt=200;
pmx = logspace(-5,100,1000);
rayl = raylpdf(pmx,1/sqrt(2));

pm_rayleigh = gather_bbechopdf_curves(N,model_dir_rayleigh,ns,glen,'rayleigh',npt,pmx);
pm_prosph = gather_bbechopdf_curves(N,model_dir_prosph,ns,glen,'prosph_0.05',npt,pmx);
pm_fish = gather_bbechopdf_curves(N,model_dir_fish,ns,glen,'fish_defaultLenDistr_defaultAngleDistr',npt,pmx);

save_folder = '/mnt/storage/broadband_echopdf_ms_figs/fig_echopdf_diff_indiv_scatterer_20141202';
if ~exist(save_folder,'dir')
    mkdir(save_folder);
end

% FIGURE: rayleigh and prolate spheroid scatterer comparison
for iN=1:length(N)
    figure;
    hh = loglog(pmx,rayl,'color',[0.5,0.5,0.5]);
    hold on
    hray = loglog(pmx,pm_rayleigh(:,iN),'k');
    hfish = loglog(pmx,pm_fish(:,iN),'b');
    hpro = loglog(pmx,pm_prosph(:,iN),'r');
    axis([1e-2 5e1 1e-6 5e1])
    title(sprintf('N=%d',N(iN)));
    xlabel('Normalized echo amplitude');
    ylabel('PDF');
    legend([hh,hray,hfish,hpro],...
           'Rayleigh distr',...
           'Rayleigh scatterer',...
           'Fish, default angle and length distr',...
           'Prolate spheroid, ar=0.05, single length, all angles');
    saveas(gcf,sprintf('%s/cmp_echopdf_fishDefaultLenDefaultAngle_prosph0.05_N_%d.fig', ...
                       save_folder,N(iN)),'fig');
    print(sprintf('%s/cmp_echopdf_fishDefaultLenDefaultAngle_prosph0.05_N_%d.png',save_folder,N(iN)),'-dpng');
    
end

