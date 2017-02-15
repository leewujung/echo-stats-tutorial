% 2016 08 09  Compare results for phasor summation for mixed assemblages

%addpath '/mnt/storage/modeling_code_current/'
%addpath '/mnt/storage/analysis_code_current/'
addpath '~/wjl/misc_matlab_code/saveSameSize'

base_path = '~/wjl/echo_stat_figs';

% Make save path
str = strsplit(mfilename('fullpath'),'/');
str = str{end};
save_path = fullfile(base_path,str);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Set param
set_all = [1,2];
ka_all = [5,8,14];
M = 20;
Ns_all = [1,10,50];
Nw = 100;
pingnum_str_all = {'1e5','1e6','1e7'};

npt = 120;  % number of points for pe kde estimation


% Set simulated data path
simu_path = 'simulation_mixed_assem_loop_20160809_set1';

% Plot and cmp
for iS=1:length(pingnum_str_all)
for iKA=1:length(ka_all)
for iNs=1:length(Ns_all)

    pingnum_str = pingnum_str_all{iS};
    pingnum = eval(pingnum_str);
    save_fname = sprintf('%s_Ns%02d_ka%02dpi_smplN%s',str,Ns_all(iNs), ...
                         ka_all(iKA),pingnum_str);

    fig = figure;
    xr = logspace(-3,log10(2000),500);  % standard
    rayl = raylpdf(xr,1/sqrt(2));
    loglog(xr,rayl,'k','linewidth',1);
    hold on
    
    for iSet=1:length(set_all)
        simu_path = sprintf('simulation_mixed_assem_loop_20160809_set%d',set_all(iSet));
        simu_file = sprintf('pnum_%s_ka%dpi_M%d_Nw%d_Ns%d.mat',...
                            pingnum_str,ka_all(iKA),M,Nw,Ns_all(iNs));
        E = load(fullfile(base_path,simu_path,simu_file));
        %[x,p_x] = findEchoDist(E.env/sqrt(mean(E.env.^2)),npt);
        [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
        loglog(x,p_x,'linewidth',1);
    end
    
    xlabel('Normalized echo amplitude','fontsize',16);
    ylabel('PDF');
    title(sprintf('ka=%dpi, Ns=%d, smplN=%s',...
                  ka_all(iKA),Ns_all(iNs),pingnum_str),'fontsize',18);
    ll = legend('Rayleigh','set 1','set 2','fontsize',16);
    set(ll,'fontsize',16)
    set(gca,'fontsize',14);
    xlim([1e-3 1e2]);
    ylim([1e-6 1e3]);

    saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
    saveSameSize_100(fig,'file',[fullfile(save_path,save_fname),'.png'],...  
                     'format','png');

end
end
end

