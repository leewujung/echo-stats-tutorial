

param.c = 1500;
param.gate_len = 0.5;

param.bp_path = '/Volumes/wjlee_apl_2/echo_stat_tutorial/echo_stat_figs/make_bpf_pool_20170417/';
param.bp_file = 'bpf_a0.211m_dtheta0.010pi_fmax1500kHz_df100Hz.mat';
param.fish_path = '~/code/echo_stat_tutorial/broadband_code_current/fish_info/';
param.fish_file = 'fish_scat_response_angle-90to90deg_len19to29cm.mat';
param.fish_len_path = '~/code/echo_stat_tutorial/broadband_code_current/fish_info/';
param.fish_len_file = 'fish_len_dist.mat';

param.scatterer.type = 'fish';
param.scatterer.nbwb = 'wb';

N = [25,250,2500];
mix_r = 1;
num_sample_str = '1e3';
num_sample = eval(num_sample_str);


save_base_path = '/Volumes/wjlee_apl_2/echo_stat_tutorial/echo_stat_figs/';
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(save_base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end
save_file = 'rayleigh';

param.save_path = save_path;
param.save_file = save_file;

for iN=1:length(N)
    bbechopdf_20170417(N(iN),mix_r,num_sample,num_sample_str,param);
end

