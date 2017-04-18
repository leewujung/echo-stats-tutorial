% Create a look up table for bpir
% 2011 12 26  Use the new, correct 2-way bpir (20111208)
% 2013 07 27  Update for using the time spacing of decimated transmit
%             signal, primary change: transpose bpir_all
% 2013 07 29  Make a pool of beampattern response in frequency domain
% 2017 04 17  Beampattern response in frequency domain for 3 deg beam
%             use ka=44.2511 from 'fig_12_pb_ka_ka_num.mat'

base_path = '/Volumes/wjlee_apl_2 1/echo_stat_tutorial/echo_stat_figs/';

% Make save path
str = strsplit(mfilename('fullpath'),'/');
str = str{end};
save_path = fullfile(base_path,str);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Calculation
fmax = 1.5e6;  % [Hz]
df = 100;  % [Hz]
dtheta = 0.001;
theta = (0:dtheta:1)*pi/2;
% a = 0.054;  % ka used in thesis--AirMarLow
a= 0.21128; % ka for 3 deg beam at 50kHz

freq_bp = 0:df:fmax;
bp = bpf_2way_fcn(theta,freq_bp,a);

fname = sprintf('bpf_a%2.3fm_dtheta%2.3fpi_fmax%dkHz_df%dHz.mat',a,dtheta, ...
                fmax/1e3,df);
save([save_path,'/',fname],'freq_bp','bp','theta','a','-MAT');



