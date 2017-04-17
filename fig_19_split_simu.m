% 2016 10 26  Simulation for interspersed echo pdf
% 2016 11 01  Simulation for split echo pdf
% 2017 04 13  Update figure legend, axis labels, and curve style

% addpath '~/wjl/misc_matlab_code/saveSameSize'
addpath '~/Dropbox/0_CODE'/MATLAB/saveSameSize/

% base_path = '~/wjl/echo_stat_figs';
base_path = '/Volumes/wjlee_apl_2/echo_stat_tutorial/echo_stat_figs/';

% Make save path
str = strsplit(mfilename('fullpath'),'/');
str = str{end};
save_path = fullfile(base_path,str);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Set param
X = load('fig_12_pb_ka_ka_num.mat');
ka = X.ka_3deg;
%ka = 2*pi;

A = 0.05;
% M = 5;
M = 20;
Nw_all = [2500];
Ns_all = [25,250,2500];%,20,50,500];
pingnum_str = '1e7';
pingnum = eval(pingnum_str);

npt = 120;  % number of points for pe kde estimation


% Set simulated data path
simu_path = 'simulation_mixture_loop_20161101_set1';

% Plot and cmp
for iNw=1:length(Nw_all)
for iNs=1:length(Ns_all)
    disp(sprintf('Ns=%d, Nw=%d',Ns_all(iNs),Nw_all(iNw)));
    save_fname = sprintf('%s_A%2.2f_M%02d_Nw%04d_Ns%04d_smpl%s',...
                         str,A,M,Nw_all(iNw),Ns_all(iNs),pingnum_str);

    fig = figure;
    xr = logspace(-3,log10(2000),500);  % standard
    rayl = raylpdf(xr,1/sqrt(2));
    hr = loglog(xr,rayl,'k','linewidth',2);
    hold on

    simu_file = sprintf('pnum_%s_ka%2.4f_A%2.2f_M%02d_Nw%04d_Ns%04d.mat',...
        pingnum_str,ka,A,M,Nw_all(iNw),Ns_all(iNs));
    E = load(fullfile(base_path,simu_path,simu_file));
    %[x,p_x] = findEchoDist(E.env/sqrt(mean(E.env.^2)),npt);
    [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
    hh = loglog(x,p_x,'r','linewidth',2);
    
%     title(sprintf('A=%2.2f ,Ns=%d, Nw=%d, smplN=%s',...
%                   A,Ns_all(iNs),Nw_all(iNw),pingnum_str),...
%           'fontsize',18);
    set(gca,'fontsize',16)
    xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
    ylabel('$p_e(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
    switch iNs
        case 1
            ll = legend(hh,'Ns = 25 (0.0937)');
        case 2
            ll = legend(hh,'Ns = 250 (0.937)');
        case 3
            ll = legend(hh,'Ns = 2500 (9.37)');
    end
    set(ll,'fontsize',22);
    xlim([1e-3 1e2]);
    ylim([1e-6 1e3]);

    saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
    saveSameSize_100(fig,'file',[fullfile(save_path,save_fname),'.png'],...  
                     'format','png');

end

end
