% 2017 04 17  Compare broadband and narrowband echo pdf generated using
%             Rayleigh scatterers

clear

addpath ~/code/echo_stat_tutorial/broadband_code_current/
addpath '~/Dropbox/0_CODE/MATLAB/saveSameSize/'

% Set params
N = [25,250,2500];
mix_r = ones(length(N),1);
num_sample_str = '5e5';
num_sample = eval(num_sample_str);

save_base_path = '/Volumes/wjlee_apl_2 1/echo_stat_tutorial/echo_stat_figs/';
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(save_base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end
save_file_pre = 'rayleigh';

param.scatterer.type = 'rayl';
param.nbwb = 'nb';

param.c = 1500;
param.gate_len = 0.5;
param.bp_path = '/Volumes/wjlee_apl_2 1/echo_stat_tutorial/echo_stat_figs/make_bpf_pool_20170417/';
% param.bp_file = 'bpf_a0.211m_dtheta0.010pi_fmax1500kHz_df100Hz.mat';
param.bp_file = 'bpf_a0.211m_dtheta0.001pi_fmax1500kHz_df100Hz.mat';
param.fish_path = '~/code/echo_stat_tutorial/broadband_code_current/fish_info/';
param.fish_file = 'fish_scat_response_angle-90to90deg_len19to29cm.mat';
param.fish_len_path = '~/code/echo_stat_tutorial/broadband_code_current/fish_info/';
param.fish_len_file = 'fish_len_dist.mat';

if 0

% Run bbechopdf code
for iN=1:length(N)
    [s, param] = bbechopdf_20170417(N(iN),mix_r,num_sample,param);
    
    if strcmp(param.nbwb,'wb')
        sfname = sprintf('%s_N_%04d_r_%02d_pnum%s_glen%2.2f_wb.mat',...
            save_file_pre,N(iN),mix_r(iN),num_sample_str,param.gate_len);
    elseif strcmp(param.nbwb,'nb')
        sfname = sprintf('%s_N_%04d_r_%02d_pnum%s_glen%2.2f_nb.mat',...
            save_file_pre,N(iN),mix_r(iN),num_sample_str,param.gate_len);
    end
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
        script_name,N(iN),mix_r(iN),num_sample_str,param.gate_len);
    saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
    saveSameSize_100(fig,'file',[fullfile(save_path,save_fname),'.png'],...
        'format','png');

end

