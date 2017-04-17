% 2012 02 26  test the effect of amplitude and phase independence on the
%             amplitude statistics
% 2012 02 27  add in theta dependence to control the pulse height and width
% 2012 02 29  change from amplitude statistics to envelope statistics

function amp_phase_theta(N,sampleN,indep,s_distr,env_opt)
% INPUT
%   N         number of scatterers
%   sampleN   sample size
%   indep     whether amplitude and phase are independent
%   s_distr   whether the scatterers are delta-distributed (0)
%                                     or Rayleigh-distributed (1)
%   env_opt   envelope statistics (1)
%             amplitude statistics (0)

%sampleN = 10;
%N = 2;
%indep = 0;
%s_distr = 0;
%env_opt = 0;

Gfolder = '/mnt/storage/ECHO_STAT/20120229_amp_phase_theta_env';
Gfile = ['gauss_ir_dtheta0.001pi_freq50kHz_fs1000kHz',...
         '_minCycleNum1_maxCycleNum60.mat'];
G = load([Gfolder,'/',Gfile]);
sigHalfL = round((size(G.win_all,2)-1)/2);
sigL = size(G.win_all,2);

minCycleNum = G.minCycleNum;
maxCycleNum = G.maxCycleNum;
freq = G.freq;  % [Hz]
fs = G.fs;   % sampling freq [Hz]
t_win = G.t_all;  % win time stamp [sec]
gateCycleNum = maxCycleNum*3;

result_dir = '/mnt/storage/ECHO_STAT/20120229_amp_phase_theta_env';
if indep==0  % amp and phase not independent
    if env_opt==1
        result_file = sprintf(['not_indep_N%d_minCycleNum%d_maxCycleNum%d',...
                            '_smplNum%d_gateCycleNum%d_env'],...
                              N,minCycleNum,maxCycleNum,sampleN, ...
                              gateCycleNum);
    else
        result_file = sprintf(['not_indep_N%d_minCycleNum%d_maxCycleNum%d',...
                            '_smplNum%d_gateCycleNum%d'],...
                              N,minCycleNum,maxCycleNum,sampleN, ...
                              gateCycleNum);
    end
else  % amp and phase independet
    if env_opt==1
        result_file = sprintf(['indep_N%d_minCycleNum%d_maxCycleNum%d',...
                            '_smplNum%d_gateCycleNum%d_env'],...
                              N,minCycleNum,maxCycleNum,sampleN,gateCycleNum);
    else
        result_file = sprintf(['indep_N%d_minCycleNum%d_maxCycleNum%d',...
                            '_smplNum%d_gateCycleNum%d'],...
                              N,minCycleNum,maxCycleNum,sampleN,gateCycleNum);
    end
end
if s_distr==0  % scatterer is delta-distributed
    result_file = [result_file,'.mat'];
else  % scatter is Rayleigh-distributed
    result_file = [result_file,'_sRayleigh.mat'];
end

gateCycleL = ceil(gateCycleNum*1/freq*fs);
t_resp = (1:gateCycleL+sigL)/fs;
smpl_pt = round(size(t_resp,2)/2); % take the mid point
s = zeros(sampleN,1);
%resp = zeros(sampleN,gateCycleL+sigL);
if indep==0  % amp and phase not indep
    parfor iS=1:sampleN 
        resp_temp = zeros(1,gateCycleL+sigL);
        theta = rand(1,N)*pi/2;
        [~,idx] = min(abs( repmat(theta',1,length(G.theta))-...
                           repmat(G.theta,length(theta),1) ), [], 2);
        time = round(rand(1,N)*gateCycleL)+sigHalfL;
        if s_distr==0
            r = ones(1,N);
        else
            r = raylrnd(ones(1,N)/sqrt(2));
        end
        for iN=1:N
            sig = r(iN)*G.y_all(idx(iN),:);
            temp = [zeros(1,time(iN)-sigHalfL),sig,...
                    zeros(1,length(resp_temp)-time(iN)-sigHalfL-1)];
            resp_temp = resp_temp+temp;
        end
        if env_opt==1
            resp_temp = abs(hilbert(resp_temp));
        end
        %resp(iS,:) = resp_temp;
        s(iS) = resp_temp(smpl_pt);
    end
else  % amp and phase indep
    parfor iS=1:sampleN
        resp_temp = zeros(1,gateCycleL+sigL);
        theta = rand(1,N)*pi/2;
        [~,idx] = min(abs( repmat(theta',1,length(G.theta))-...
                           repmat(G.theta,length(theta),1) ), [], 2);
        time = round(rand(1,N)*gateCycleL)+sigHalfL;
        ph = rand(1,N)*2*pi; % random phase
        if s_distr==0
            r = ones(1,N);
        else
            r = raylrnd(ones(1,N)/sqrt(2));
        end
        for iN=1:N
            sig = r(iN)*cos(2*pi*freq*t_win+ph(iN)).*G.win_all(idx(iN),:);
            % sig = ph(iN)*r(iN)*G.win_all(idx(iN),:);
            temp = [zeros(1,time(iN)-sigHalfL),sig,...
                    zeros(1,length(resp_temp)-time(iN)-sigHalfL-1)];
            resp_temp = resp_temp+temp;
        end
        if env_opt==1
            resp_temp = abs(hilbert(resp_temp));
        end
        %resp(iS,:) = resp_temp;
        s(iS) = resp_temp(smpl_pt);
    end    
end
save([result_dir,'/',result_file],...
     's','fs','sampleN','N','minCycleNum','maxCycleNum',...
     'gateCycleNum','indep','s_distr');
