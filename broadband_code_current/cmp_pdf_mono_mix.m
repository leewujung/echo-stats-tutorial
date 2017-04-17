% 2012 05 14  Compare two sets of modeling result
% 2012 08 02  Update so that more than 2 sets of files can be compared
% 2012 11 12  Compare monotype and mixed assemblage pdf model

clear t

addpath '/mnt/storage/modeling_code_current/'
addpath '/mnt/storage/analysis_code_current/'

save_dir = '/mnt/storage/ECHO_STAT/20121112_mono_mix_cmp';
save_opt = 1;

sampleN = 5e4;
Nscat_all = [10,50,300];

folder{1} = '/mnt/storage/ECHO_STAT/20121108_bbfig_models';
folder{2} = '/mnt/storage/ECHO_STAT/20121108_bbfig_models';

fname_pre = 'mono_mix_cmp';

% Filename prefix
t{1} = 'tx_sys';
t{2} = 'mixed_N_1000';

% Filename postfix
te{1} = '_freqDepBP.mat';
te{2} = '_freqDepBP.mat';


r = 30;  % ratio of strong to weak scatterer
t_len = 'Nw=1000, Ns=';

tc = {'b','r','k','g','m','c'};
ts = {'o','x','+','^','v','*'};

tail = ',''markersize'',6,''linewidth'',1);';
xr = logspace(-3,log10(2000),500);  % standard
rayl = raylpdf(xr,1/sqrt(2));
npt = 1e2;
for N = Nscat_all
    ff{1} = [t{1},'_N',num2str(N),'_sampleN',num2str(sampleN),...
                 '_gateLdist0.5',te{1}];
    ff{2} = [t{2},'_',num2str(N),'_r_1_',num2str(r),'_sampleN',num2str(sampleN),...
                 '_gateLen0.5',te{2}];
    for iF=1:2
    str{iF} = [folder{iF},'/',ff{iF}];
    s{iF} = load(str{iF});
    [dens(iF,:),x(iF,:),bw(iF)] =...
        findEchoDist_kde(s{iF}.s/sqrt(mean(s{iF}.s.^2)),npt);
    end

    fig = figure;
    loglog(xr,rayl,'k');
    hold on
    loglog(x(1,:),dens(1,:),tc{1});
    loglog(x(2,:),dens(2,:),tc{2});
    xlim([1e-3 50]);
    ylim([1e-5 1e2]);
    legend('Rayleigh','Monotype','Mixed assemblage','location','southwest');
    set(gca,'fontsize',16);
    tt = [t_len,num2str(N),', r=',num2str(r)];
    title(tt,'fontsize',20);
    if save_opt==1
        save_fname = [t{2},'_',num2str(N),'_r',num2str(r),'_mono_mix_cmp'];
        saveas(gca,[save_dir,'/',save_fname,'.fig'],'fig');
        saveas(gca,[save_dir,'/',save_fname,'.png'],'png');
    end

end

