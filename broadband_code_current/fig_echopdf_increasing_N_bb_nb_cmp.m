% 2014 06 26  Plot figure to show the non-Rayleigh to Rayleigh
%             distribution with increasing N
% 2014 09 01  Add in narrowband echo pdf with increasing N

% Get bbechopdf model for plotting the pool
%N=[10,20,50,100,300,1000]';
N=[10,50,300]';
ns = 5e4;
glen = 0.5;
model_dir_prosph   = '/mnt/storage/broadband_echopdf_models/prosph_0.5';
model_dir_rayleigh = '/mnt/storage/broadband_echopdf_models/rayleigh_scatterers';
model_dir_fish     = '/mnt/storage/broadband_echopdf_models/fish_defaultLenDistr_defaultAngleDistr';
model_dir_nb       = '/mnt/storage/broadband_echopdf_models/narrowband_rayleigh';

npt=200;
pmx = logspace(-5,100,1000);
rayl = raylpdf(pmx,1/sqrt(2));

pm_rayleigh = gather_bbechopdf_curves(N,model_dir_rayleigh,ns,glen,'rayleigh',npt,pmx);
pm_prosph   = gather_bbechopdf_curves(N,model_dir_prosph,ns,glen,'prosph_0.5',npt,pmx);
pm_fish     = gather_bbechopdf_curves(N,model_dir_fish,ns,glen,'fish_defaultLenDistr_defaultAngleDistr',npt,pmx);
pm_nb       = gather_bbechopdf_curves_nb(N,model_dir_nb,ns,glen,'rayleigh',npt,pmx);

save_folder = '/mnt/storage/broadband_echopdf_ms_figs/fig_echopdf_increasing_N_20140901';
if ~exist(save_folder,'dir')
    mkdir(save_folder);
end


%% FIGURE: individual comparison between broadband and narrowband
%% echo pdf, both use Rayleigh scatterers

for iN=1:length(N)

figure;
loglog(pmx,rayl,'color',[0.5,0.5,0.5]);
hold on
loglog(pmx,pm_nb(:,iN),'b');
loglog(pmx,pm_rayleigh(:,iN),'r');
axis([1e-2 5e1 5e-5 5e1])
title(['Rayleigh scatterer, bb vs nb, N=',num2str(N(iN))]);
xlabel('Normalized echo amplitude');
ylabel('PDF');
legend('Rayleigh','narrowband','broadband','location','southwest');
saveas(gcf,sprintf('%s/bb_vs_nb_N%d_rayleigh_scatterer.png', ...
                   save_folder,N(iN)),'png');
saveas(gcf,sprintf('%s/bb_vs_nb_N%d_rayleigh_scatterer.fig', ...
                   save_folder,N(iN)),'fig');

end


