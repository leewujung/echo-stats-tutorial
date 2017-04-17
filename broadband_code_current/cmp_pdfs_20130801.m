% 2012 05 14  Compare two sets of modeling result
% 2012 08 02  Update so that more than 2 sets of files can be compared

clear t

addpath '/mnt/storage/modeling_code_current/'
addpath '/mnt/storage/analysis_code_current/'

%save_dir = '/mnt/storage/ECHO_STAT/20130801_bb_pdfmodel_cmp_rayleigh_fish';
%save_dir = '/mnt/storage/ECHO_STAT/20130801_bb_pdfmodel_rayleigh_cmp_old_new';
Nscat_all = [10:10:100,200:100:1000];

folder{1} = '/mnt/storage/ECHO_STAT/20130801_bb_pdfmodel_rayleigh';
folder{2} = '/mnt/storage/ECHO_STAT/20121108_pdfmodel_sample_only';
%folder{2} = '/mnt/storage/ECHO_STAT/20130801_bb_pdfmodel_fish';

fname_pre = 'rayleigh_fish_indiv';

% Sample number
sn{1} = 5e4;
sn{2} = 5e4;

% Filename prefix
t{1} = 'tx_sys';
t{2} = 'tx_sys';

% Filename postfix
te{1} = '_freqDepBP.mat';
te{2} = '_freqDepBP.mat';

% Legend
tt{1} = 'Rayleigh indiv';
tt{2} = 'fish indiv';

tc = {'r','b','k','g','m','c'};
ts = {'o','x','+','^','v','*'};

tail = ',''markersize'',6,''linewidth'',1);';
M{1} = 'Rayleigh';
for N = Nscat_all
    fig = figure;
    xr = logspace(-3,log10(2000),500);  % standard
    rayl = raylpdf(xr,1/sqrt(2));
    loglog(xr,rayl,'k');
    hold on
    for iF = 1:length(t)
        if iF==1
            fname = [t{iF},'_N',num2str(N),'_sampleN',num2str(sn{iF}),...
                     '_gateLen0.5',te{iF}];
        elseif iF==2
            fname = [t{iF},'_N',num2str(N),'_sampleN',num2str(sn{iF}),...
                     '_gateLdist0.5',te{iF}];
        end
        str = [folder{iF},'/',fname];
        eval(['s',num2str(iF),'=load(''',str,''',''s''',');']);
        eval(['[p_s',num2str(iF)',',x] = findEchoDist_kde(s',num2str(iF),'.s./sqrt(mean(s',num2str(iF),'.s.^2)),100);']);
        rayl = raylpdf(x,1/sqrt(2));
        figure(fig);
        %eval(['loglog(x,p_s',num2str(iF),',''',tc{iF},ts{iF},'''',tail]);
        eval(['loglog(x,p_s',num2str(iF),',''',tc{iF},'''',tail]);
        %ttemp = t{iF};
        %ttemp(ttemp=='_') = ' ';
        M{iF+1} = tt{iF};
    end
    xlim([1e-3 1e2]);
    ylim([1e-7 1e2]);
    grid on
    legend(M,'location','southwest');
    set(gca,'fontsize',16);
    title(['N=',num2str(N)],'fontsize',20);
    save_fname = [fname_pre,'_N',num2str(N)];
    saveas(gca,[save_dir,'/',save_fname,'.fig'],'fig');
    saveas(gca,[save_dir,'/',save_fname,'.png'],'png');
    close
end

