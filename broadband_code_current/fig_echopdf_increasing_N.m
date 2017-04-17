% 2014 06 26  Plot figure to show the non-Rayleigh to Rayleigh
%             distribution with increasing N
% 2014 09 01  Add in narrowband echo pdf with increasing N
% 2014 11 23  Change the prolate spheroid scatterer echopdf models
%             to those generated using aspect_ratio=0.05 and
%             default length and angle of orientation distribution

% Get bbechopdf model for plotting the pool
%N=[10,20,50,100,300,1000]';
N=[10,20,50,100,200,500,800]';
ns = 5e4;
glen = 0.5;
model_dir_prosph   = '/mnt/storage/broadband_echopdf_models/prosph_0.05_defaultAngleDistr_defaultLenDistr';
model_dir_rayleigh = '/mnt/storage/broadband_echopdf_models/rayleigh_scatterers';
model_dir_fish     = '/mnt/storage/broadband_echopdf_models/fish_defaultLenDistr_defaultAngleDistr';
model_dir_nb       = '/mnt/storage/broadband_echopdf_models/narrowband_rayleigh';

npt=200;
pmx = logspace(-5,100,1000);
rayl = raylpdf(pmx,1/sqrt(2));

pm_rayleigh = gather_bbechopdf_curves(N,model_dir_rayleigh,ns,glen,'rayleigh',npt,pmx);
pm_prosph   = gather_bbechopdf_curves(N,model_dir_prosph,ns,glen,'prosph_0.05',npt,pmx);
pm_fish     = gather_bbechopdf_curves(N,model_dir_fish,ns,glen,'fish_defaultLenDistr_defaultAngleDistr',npt,pmx);
pm_nb       = gather_bbechopdf_curves_nb(N,model_dir_nb,ns,glen,'rayleigh',npt,pmx);

save_folder = '/mnt/storage/broadband_echopdf_ms_figs/fig_echopdf_increasing_N_20141123';
if ~exist(save_folder,'dir')
    mkdir(save_folder);
end

% FIGURE: increasing N --> echopdf toward Rayleigh
% narrowband Rayleigh scatterer
figure;
loglog(pmx,rayl,'color',[0.5,0.5,0.5]);
hold on
loglog(pmx,pm_nb,'b');
axis([1e-2 5e1 5e-5 5e1])
title('Rayleigh scatterers, narrowband');
xlabel('Normalized echo amplitude');
ylabel('PDF');
legend({'Rayleigh',num2str([N])});
saveas(gcf,sprintf('%s/narrowband_rayleigh_scatterer_increasingN.png', ...
                   save_folder),'png');
saveas(gcf,sprintf('%s/narrowband_rayleigh_scatterer_increasingN.fig', ...
                   save_folder),'fig');

% broadband Rayleigh scatterer
figure;
loglog(pmx,rayl,'color',[0.5,0.5,0.5]);
hold on
loglog(pmx,pm_rayleigh,'b');
axis([1e-2 5e1 5e-5 5e1])
title('Rayleigh scatterers, broadband');
xlabel('Normalized echo amplitude');
ylabel('PDF');
legend({'Rayleigh',num2str([N])});
saveas(gcf,sprintf('%s/rayleigh_scatterer_increasingN.png', ...
                   save_folder),'png');
saveas(gcf,sprintf('%s/rayleigh_scatterer_increasingN.fig', ...
                   save_folder),'fig');

% Broadband prolate spheroid scatterer
figure;
loglog(pmx,rayl,'color',[0.5,0.5,0.5]);
hold on
loglog(pmx,pm_prosph,'b');
axis([1e-2 5e1 5e-5 5e1])
title('Prolate spheroid scatterers, broadband');
xlabel('Normalized echo amplitude');
ylabel('PDF');
legend({'Rayleigh',num2str([N])});
saveas(gcf,sprintf('%s/prolate_spheroid_scatterer_increasingN.png', ...
                   save_folder),'png');
saveas(gcf,sprintf('%s/prolate_spheroid_scatterer_increasingN.fig', ...
                   save_folder),'fig');

% broadband fish scatterer
figure;
loglog(pmx,rayl,'color',[0.5,0.5,0.5]);
hold on
loglog(pmx,pm_fish,'b');
axis([1e-2 5e1 5e-5 5e1])
title('Fish scatterers, broadband');
xlabel('Normalized echo amplitude');
ylabel('PDF');
legend({'Rayleigh',num2str([N])});
saveas(gcf,sprintf('%s/fish_scatterer_increasingN.png', ...
                   save_folder),'png');
saveas(gcf,sprintf('%s/fish_scatterer_increasingN.fig', ...
                   save_folder),'fig');


%% FIGURE: individual comparison between broadband and narrowband
%% echo pdf, both use Rayleigh scatterers --> move to another file



