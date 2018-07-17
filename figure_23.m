% Code to generate Figure 23 of the echo statistics tutorial.
%
% This code plots comparisons between PDFs of the echo magnitudes
% associated with long narrowband and short-pulsed broadband signals.
% The N identical Rayleigh scatterers are randomly distributed in the
% sensor beam.
% 3D distribution of scatterers.
%
% Package 'broadband-echo-stats' is used to generate the ensemble of
% samples for both narrowband and broadband echo pdfs for this comparison
% figure. This package is located at: https://github.com/leewujung/broadband-echo-stats
% If you have cloned this code from GitHub this package is included as a
% submodule. To make sure you have the code, run
%     git submodule init
% and then
%     git submodule update
%
% Author: Wu-Jung Lee | leewujung@gmail.com | APL-UW


% 2017 04 17  Compare broadband and narrowband echo pdf generated using
%             Rayleigh scatterers

clear
addpath './util_fcn'
addpath './broadband-echo-stats'
base_path = './figs';
bpf_path = './figs/make_bpf_pool';
fish_info_path = './broadband-echo-stats/fish_info';

% Make save path
str = strsplit(mfilename('fullpath'),'/');
str = str{end};
save_path = fullfile(base_path,str);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Test if already have broadband beampattern response lookup table
% if not, run code to generate one
if ~exist(bpf_path,'dir')
    disp('Do not have broadband beampattern response lookup table.');
    disp('Generating one...');
    make_bpf_pool;
end

% Set overall params
N = [25,250,2500];
mix_r = ones(length(N),1);
num_sample_str = '1e6';
num_sample = eval(num_sample_str);
save_file_pre = 'rayleigh';

% Set nbwb code params
param.scatterer.type = 'rayl';
param.c = 1500;
param.gate_len = 0.5;
param.bp_path = bpf_path;
param.bp_file = 'bpf_a0.211m_dtheta0.001pi_fmax1500kHz_df100Hz.mat';
param.fish_path = fish_info_path;
param.fish_file = 'fish_scat_response_angle-90to90deg_len19to29cm.mat';
param.fish_len_path = fish_info_path;
param.fish_len_file = 'fish_len_dist.mat';


% Set operation
mc_opt = 1;  % 0 - do not re-generate realizations
             % 1 - re-generate all realizations

% Monte Carlo simulation
if mc_opt
    % Run bbechopdf code
    for iN=1:length(N)
        disp('Generating broadband samples');
        param.nbwb = 'wb';
        [s, param] = bbechopdf_20170417(N(iN),mix_r,num_sample,param);
        sfname = sprintf('%s_N_%04d_r_%02d_pnum%s_glen%2.2f_wb.mat',...
                save_file_pre,N(iN),mix_r(iN),num_sample_str,param.gate_len);
        save([save_path,'/',sfname],'param','s');

        disp('Generating narrowband samples');
        param.nbwb = 'nb';
        [s, param] = bbechopdf_20170417(N(iN),mix_r,num_sample,param);
        sfname = sprintf('%s_N_%04d_r_%02d_pnum%s_glen%2.2f_nb.mat',...
                save_file_pre,N(iN),mix_r(iN),num_sample_str,param.gate_len);
        save([save_path,'/',sfname],'param','s');
    end
end



% Plot
npt = 120;
for iN=1:length(N)
    fig = figure;
    xr = logspace(-3,log10(2000),500);
    rayl = raylpdf(xr,1/sqrt(2));
    loglog(xr,rayl,'k','linewidth',2);
    hold on

    % Load narrowband file
    simu_file = sprintf('%s_N_%04d_r_%02d_pnum%s_glen%2.2f_nb.mat',...
        'rayleigh',N(iN),1,num_sample_str,param.gate_len);
    E = load(fullfile(save_path,simu_file));
    [p_x,x] = findEchoDist_kde(E.s/sqrt(mean(E.s.^2)),npt);
%     [x,p_x] = findEchoDist(E.s/sqrt(mean(E.s.^2)),npt);
    loglog(x,p_x,'b-','linewidth',2);

    % Load broadband file
    simu_file = sprintf('%s_N_%04d_r_%02d_pnum%s_glen%2.2f_wb.mat',...
        'rayleigh',N(iN),1,num_sample_str,param.gate_len);
    E = load(fullfile(save_path,simu_file));
    [p_x,x] = findEchoDist_kde(E.s/sqrt(mean(E.s.^2)),npt);
%     [x,p_x] = findEchoDist(E.s/sqrt(mean(E.s.^2)),npt);
    loglog(x,p_x,'r-','linewidth',2);

    % misc
    title(sprintf('N=%d, a=0.211m, smplN=%s',N(iN),num_sample_str),...
        'fontsize',18);

    ll = legend('Rayleigh','Narrowband','Broadband',...
        'location','southwest');
    set(ll,'fontsize',18);
    set(gca,'fontsize',16)
    xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
    ylabel('$p_e(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
    xlim([1e-3 1e2]);
    ylim([1e-6 1e3]);

    % Save figure
    save_fname = sprintf('%s_N_%04d_r_%02d_pnum%s_glen%2.2f',...
        str,N(iN),mix_r(iN),num_sample_str,param.gate_len);
    saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
    saveSameSize(fig,'file',[fullfile(save_path,save_fname),'.png'],...
        'format','png');

end
