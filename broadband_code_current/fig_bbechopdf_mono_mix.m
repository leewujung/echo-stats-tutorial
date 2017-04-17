% 2014 06 30  Plot broadband echo pdfs generated assuming monotype
%             and mixed assemblages

addpath /mnt/storage/analysis_code_current

MONO_DIR = '/mnt/storage/broadband_echopdf_models/rayleigh_scatterers';
MIX_DIR = '/mnt/storage/broadband_echopdf_models/mixed_assemblage';
SAVE_DIR = '/mnt/storage/broadband_echopdf_ms_figs/fig_echopdf_mono_mix_cmp';
if ~exist(SAVE_DIR,'dir')
    mkdir(SAVE_DIR);
end

r = [5,10,30];
Ndom = [10,50,300];
ns = 5e4;
glen = 0.5;
tail = sprintf('_sampleN%d_gateLen%2.1f_freqDepBP.mat',...
               ns,glen);
npt = 200;

pmx = logspace(-5,100,3000);
rayl = raylpdf(pmx,1/sqrt(2));

for iR=1:length(r)
for iN=1:length(Ndom)
    fname_mono = ['rayleigh','_N_',num2str(Ndom(iN)),'_r_1',tail];
    A_mono = load([MONO_DIR,'/',fname_mono]);
    [pm_mono,pmx_mono] = findEchoDist_kde(A_mono.s/sqrt(mean(A_mono.s.^2)),npt);
    pm_mono(isnan(pm_mono)==1) = [];
    pmx_mono(isnan(pm_mono)==1) = [];    
    pm_mono_final = interp1(pmx_mono,pm_mono,pmx);

    fname_mix = ['rayleigh','_N_1000_',num2str(Ndom(iN)),'_r_1_',num2str(r(iR)),tail];
    A_mix = load([MIX_DIR,'/',fname_mix]);
    [pm_mix,pmx_mix] = findEchoDist_kde(A_mix.s/sqrt(mean(A_mix.s.^2)),npt);
    pm_mix(isnan(pm_mix)==1) = [];
    pmx_mix(isnan(pm_mix)==1) = [];    
    pm_mix_final = interp1(pmx_mix,pm_mix,pmx);
    
    figure;
    loglog(pmx,rayl,'color',[0.5,0.5,0.5]);
    hold on
    loglog(pmx,pm_mono_final,'b');
    loglog(pmx,pm_mix_final,'r');
    axis([1e-2 5e1 1e-6 5e1])
    title(sprintf('N=%d, r=%d',Ndom(iN),r(iR)));
    xlabel('Normalized echo amplitude');
    ylabel('PDF');
    legend('Rayleigh distr','Mono','Mix');
    saveas(gcf,sprintf('%s/rayleigh_scatterers_3000pt_mix_r_%d_Ndom_%d.fig',...
                       SAVE_DIR,r(iR),Ndom(iN)),'fig');
    print(sprintf('%s/rayleigh_scatterers_3000pt_mix_r_%d_Ndom_%d.png',...
                  SAVE_DIR,r(iR),Ndom(iN)),'-dpng');

end
end



