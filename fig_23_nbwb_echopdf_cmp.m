% 2017 04 17  Compare broadband and narrowband echo pdf generated using
%             Rayleigh scatterers

clear

addpath ~/code/echo_stat_tutorial/broadband_code_current/

% Set params
N = [25,250];
mix_r = [1,1];
num_sample_str = '1e3';
num_sample = eval(num_sample_str);

save_base_path = '/Volumes/wjlee_apl_2/echo_stat_tutorial/echo_stat_figs/';
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(save_base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end
save_file_pre = 'rayleigh';

param.scatterer.type = 'rayl';
param.scatterer.nbwb = 'wb';

param.c = 1500;
param.gate_len = 0.5;
param.bp_path = '/Volumes/wjlee_apl_2/echo_stat_tutorial/echo_stat_figs/make_bpf_pool_20170417/';
param.bp_file = 'bpf_a0.211m_dtheta0.010pi_fmax1500kHz_df100Hz.mat';
param.fish_path = '~/code/echo_stat_tutorial/broadband_code_current/fish_info/';
param.fish_file = 'fish_scat_response_angle-90to90deg_len19to29cm.mat';
param.fish_len_path = '~/code/echo_stat_tutorial/broadband_code_current/fish_info/';
param.fish_len_file = 'fish_len_dist.mat';

% Run bbechopdf code
for iN=1:length(N)
    [s, param] = bbechopdf_20170417(N(iN),mix_r,num_sample,param);
    
    sfname = sprintf('%s_N_%04d_r_%02d_pnum%s_glen%2.2f_freqDepBP.mat',...
        save_file_pre,N(iN),mix_r(iN),num_sample_str,param.gate_len);
    save([save_path,'/',sfname],'param','s');

end


% Plot
npt = 200;
for iN=1:length(N)
    fig = figure;
    xr = logspace(-3,log10(2000),500);
    rayl = raylpdf(xr,1/sqrt(2));
    loglog(xr,rayl,'k','linewidth',2);
    hold on
    
    % Load broadband file
    simu_file = sprintf('%s_N_%04d_r_%02d_pnum%s_glen%2.2f_freqDepBP.mat',...
        'rayleigh',N(iN),1,num_sample_str,param.gate_len);
    E = load(fullfile(save_path,simu_file));
    [p_x,x] = findEchoDist_kde(E.s/sqrt(mean(E.s.^2)),npt);
    loglog(x,p_x,'r-','linewidth',2);
    
    title(sprintf('N=%d, a=0.211m, smplN=%s',N(iN),num_sample_str),...
        'fontsize',18);

    ll = legend('Rayleigh','Broadband',...
        'location','southwest');
    set(ll,'fontsize',18);
    set(gca,'fontsize',16)
    xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
    ylabel('$p_e(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
    xlim([1e-3 1e2]);
    ylim([1e-6 1e3]);

end

