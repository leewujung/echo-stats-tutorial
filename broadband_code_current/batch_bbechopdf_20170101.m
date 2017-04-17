

N = [25,250,2500];
mix_r = 1;
ns = 1e3;
gate_len = 0.5;
save_opt = 1;
tx = 1;

save_base_path = '~/Desktop/echo_stat_figs';

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(save_base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

sdir = save_path;

for iN=1:length(N)
    % bbechopdf(N,mix_r,ns,gate_len,tx,sdir,fname,...
    %           bpa,indiv,[len_bin,len_dist],[angle_mean,angle_std],nb_freq);
    
    % Rayleigh
    bbechopdf(N(iN),mix_r,ns,gate_len,tx,sdir,'rayleigh',[],0);
end

