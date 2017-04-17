% Create a look up table for bpir
% 2011 12 26  Use the new, correct 2-way bpir (20111208)
% 2013 07 27  Update for using the time spacing of decimated transmit
%             signal, also transpose bpir_all

%% Get transmit signal
tx_opt = 3;
[y,t_y] = gen_tx(tx_opt);  % signal already decimated

fmax = 1.5e6;  % [Hz]
dtheta = 0.001;
theta = (dtheta:dtheta:1)*pi/2;
a = 0.054;

[t90,bt90] = bpir_2way_fcn(pi/2,fmax,a);
[p,q] = rat(diff(t90(1:2))/diff(t_y(1:2)));
bt90_y = resample(bt90,p,q);
t90_y = (0:length(bt90_y)-1)*diff(t90(1:2))*q/p;
maxHalfL = round((length(bt90_y)+1)/2);  % half of the longest bpir

bpir = zeros(length(bt90_y),length(theta));  % all bpir [length x theta]

folder = '/mnt/storage/broadband_code_current/bpir_bpf_pool';
parfor iT=1:length(theta)
    disp(['iT=',num2str(iT)]);
    [t,bt] = bpir_2way_fcn(theta(iT),fmax,a);
    bt_y = resample(bt,p,q); % resample to match the transmit signal
    tHalfL = round((length(bt_y)+1)/2);
    if mod(length(bt_y),2)~=0
        bt_temp = [zeros(1,maxHalfL-tHalfL),bt_y,...
                   zeros(1,maxHalfL-tHalfL)];
    else
        bt_temp = [zeros(1,maxHalfL-tHalfL),bt_y,...
                   zeros(1,maxHalfL-tHalfL+1)];
    end
    bpir(:,iT) = bt_temp;
end
t_bpir = t90;

fname = sprintf('bpir_a%2.3fm_dtheta%2.3fpi_fmax%dkHz_matchTX.mat',a,dtheta, ...
                fmax/1e3);
save([folder,'/',fname],'t_bpir','bpir','theta','a',...
     'fmax','dtheta','y','t_y','-MAT');



