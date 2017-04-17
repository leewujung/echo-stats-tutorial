% 2012 05 14  Compare two sets of modeling result
% 2012 08 02  Update so that more than 2 sets of files can be compared

clear t

addpath '/mnt/storage/modeling_code_current/'
addpath '/mnt/storage/analysis_code_current/'

save_dir = '/mnt/storage/ECHO_STAT/20120814_bb_1e4smpl_2_cmp';
sampleN = 1e5;
Nscat_all = [10];

%folder{1} = '/mnt/storage/ECHO_STAT/20120731_bp_results';
%folder{2} = '/mnt/storage/ECHO_STAT/20120731_bp_results';
folder{1} = '/mnt/storage/ECHO_STAT/20120814_bb_1e4smpl_2';
folder{2} = '/mnt/storage/ECHO_STAT/20120814_bb_1e4smpl_2';
%folder{3} = '/mnt/storage/ECHO_STAT/20120731_bp_results';
%folder{4} = '/mnt/storage/ECHO_STAT/20120731_bp_results';
%folder{5} = '/mnt/storage/ECHO_STAT/20120731_bp_results';
%folder{6} = '/mnt/storage/ECHO_STAT/20120731_bp_results';

fname_pre = 'sqchrip_txsys';

% Filename prefix
%t{1} = 'middle_hann_taper_narrow';
%t{2} = 'middle_hann_taper';
%t{3} = 'middle_hann_taper_wide';
%t{4} = 'square_chirp_no_taper';
%t{5} = 'square_chirp_no_taper';
%t{6} = 'square_chirp_no_taper';
%t{1} = 'lower_hann_taper';
%t{2} = 'middle_hann_taper';
%t{3} = 'upper_hann_taper';
%t{1} = 'no_taper';
%t{2} = 'actual_w_sys_no_taper';
%t{2} = 'hann_hf_taper';
%t{3} = 'hann_lf_taper';
%t{1} = 'sqchirp_narrow1_hann';
%t{2} = 'sqchirp_wide2_hann';
t{1} = 'tx_nosys';
t{2} = 'tx_sys';


% Filename postfix
%t{3} = 'no_taper';
%t{4} = 'actual_w_sys_no_taper';
te{1} = '_freqDepBP.mat';
%te{1} = '_nobp.mat';
te{2} = '_freqDepBP.mat';
%te{2} = '_fixfreq30kHz.mat';
%te{3} = '_fixfreq40kHz.mat';
%te{4} = '_fixfreq50kHz.mat';
%te{5} = '_fixfreq60kHz.mat';
%te{6} = '_fixfreq70kHz.mat';

% Legend
%tt{1} = 'sq no taper, freq-dep';
%tt{2} = 'sq no taper, 30kHz';
%tt{3} = 'sq no taper, 40kHz';
%tt{4} = 'sq no taper, 50kHz';
%tt{5} = 'sq no taper, 60kHz';
%tt{6} = 'sq no taper, 70kHz';
%tt{1} = 'middle hann taper narrow';
%tt{2} = 'middle hann taper';
%tt{3} = 'middle hann taper wide';
%tt{1} = 'lower hann win';
%tt{2} = 'middle hann win';
%tt{3} = 'upper hann win';
%tt{1} = 'narrow Hann';
%tt{2} = 'wide Hann';
%tt{1} = 'without bp';
%tt{2} = 'with bp';
tt{1} = 'tx no sys';
tt{2} = 'tx with sys';

tc = {'b','r','k','g','m','c'};
ts = {'x','o','+','^','v','*'};

tail = ',''markersize'',6,''linewidth'',1);';
M{1} = 'Rayleigh';
for N = Nscat_all
    fig = figure;
    xr = logspace(-3,log10(2000),500);  % standard
    rayl = raylpdf(xr,1/sqrt(2));
    loglog(xr,rayl,'k');
    hold on
    for iF = 1:length(t)
        fname = [t{iF},'_N',num2str(N),'_sampleN',num2str(sampleN),...
                 '_gateLdist0.5',te{iF}];
        str = [folder{iF},'/',fname];
        eval(['s',num2str(iF),'=load(''',str,''',''s''',');']);
        eval(['[x,p_s',num2str(iF)','] = findEchoDist(s',num2str(iF),'.s,''default'',0,90);']);

        npt = 150;
        xf = logspace(-3,log10(2000),npt);
        eval(['data = s',num2str(iF),'.s;']);
        eval(['f_s',num2str(iF),'=ksdensity(data/sqrt(mean(data.^2)),xf);']);
        
        rayl = raylpdf(x,1/sqrt(2));
        figure(fig);
        %eval(['loglog(x,p_s',num2str(iF),',''',tc{iF},ts{iF},'''',tail]);
        eval(['loglog(xf,f_s',num2str(iF),',''',tc{iF},'''',tail]);
        %ttemp = t{iF};
        %ttemp(ttemp=='_') = ' ';
        M{iF+1} = tt{iF};
    end
    xlim([1e-3 1e2]);
    ylim([1e-5 1e2]);
    grid on
    legend(M,'location','southwest');
    set(gca,'fontsize',16);
    title(['N=',num2str(N)],'fontsize',20);
    save_fname = [fname_pre,'_N',num2str(N),'_line'];
    saveas(gca,[save_dir,'/',save_fname,'.fig'],'fig');
    saveas(gca,[save_dir,'/',save_fname,'.png'],'png');
    close
end

