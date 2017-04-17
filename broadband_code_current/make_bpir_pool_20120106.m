% Creat a look up table for bpir
% 2011 12 26  Use the new, correct 2-way bpir (20111208)

fmax = 1e7;
dtheta = 0.001
theta = (dtheta:dtheta:1)*pi/2;
a = 0.054;

[t90,bt90] = bpir_2way_fcn(pi/2,fmax,a);
maxHalfL = (length(t90)+1)/2;  % half of the longest bpir

bpir_all = zeros(length(theta),length(t90));  % all bpir [length x theta]

folder = '/mnt/storage/ECHO_STAT/20120223_bp_results';
parfor iT=1:length(theta)
    disp(['iT=',num2str(iT)]);
    [t,bt] = bpir_2way_fcn(theta(iT),fmax,a);
    tHalfL = (length(t)+1)/2;
    bt_temp = [zeros(1,maxHalfL-tHalfL),bt,zeros(1,maxHalfL-tHalfL)];
    bpir_all(iT,:) = bt_temp;
end

t_all = t90;

fname = sprintf('bpir_a%2.3fm_dtheta%2.3fpi_fmax%dkHz.mat',a,dtheta, ...
                fmax/1e3);
save([folder,'/',fname],'t_all','bpir_all','theta','a',...
     'fmax','dtheta','-MAT');



