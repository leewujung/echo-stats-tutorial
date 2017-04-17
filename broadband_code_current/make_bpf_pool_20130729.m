% Create a look up table for bpir
% 2011 12 26  Use the new, correct 2-way bpir (20111208)
% 2013 07 27  Update for using the time spacing of decimated transmit
%             signal, primary change: transpose bpir_all
% 2013 07 29  Make a pool of beampattern response in frequency domain

fmax = 1.5e6;  % [Hz]
df = 100;  % [Hz]
dtheta = 0.01;
theta = (dtheta:dtheta:1)*pi/2;
a = 0.054;

freq_bp = 0:df:fmax;
bp = bpf_2way_fcn(theta,freq_bp,a);

folder = '/mnt/storage/broadband_code_current/bpir_bpf_pool';
fname = sprintf('bpf_a%2.3fm_dtheta%2.3fpi_fmax%dkHz_df%dHz.mat',a,dtheta, ...
                fmax/1e3,df);
save([folder,'/',fname],'freq_bp','bp','theta','a','-MAT');



