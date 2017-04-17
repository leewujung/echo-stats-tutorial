% 2012 02 26  test the effect of amplitude and phase independence on the
%             amplitude statistics
% 2012 02 29  change from amplitude statistics to envelope statistics

function amp_phase_indep_test(N,cycleNum,sampleN,gateCycleNum,indep,s_distr)
% INPUT
%   N              number of scatterers
%   cycleNum       how many cycles are there in one pulse
%   sampleN        sample size
%   gateCycleNum   number of cycles in the gate
%   indep          whether amplitude and phase are independent
%   s_distr        whether the scatterers are delta-distributed (0)
%                                       or Rayleigh-distributed (1)

%cycleNum = 10;
%gateCycleNum = 50;
%sampleN = 1e5;
%N = 10;
%indep = 0;

freq = 50e3;  % [Hz]
fs = 5e6;   % sampling freq [Hz]

result_dir = '/mnt/storage/ECHO_STAT/20120229_amp_phase_indep_env';
if indep==0
    result_file = sprintf('not_indep_N%d_cycleNum%d_smplNum%d_gateCycleNum%d',...
                          N,cycleNum,sampleN,gateCycleNum);
else
    result_file = sprintf('indep_N%d_cycleNum%d_smplNum%d_gateCycleNum%d',...
                          N,cycleNum,sampleN,gateCycleNum);
end
if s_distr==0
    result_file = [result_file,'.mat'];
else
    result_file = [result_file,'_sRayleigh.mat'];
end


gateCycleL = ceil(gateCycleNum*1/freq*fs);

t_y = ((1:ceil(1/freq*cycleNum*fs))-1)/fs; % time stamp for signal
y = sin(2*pi*freq*t_y);

winL = length(y);
win = gausswin(winL);
win = [win(1:floor(winL/2))',ones(1,length(y)-winL),...
       win(floor(winL/2)+1:end)'];
if indep==0  % amp and phase not indep
    y = y.*win;
else  % amp and phase indep
    y = win;
end

t_resp = (1:gateCycleL+length(y))/fs;
%resp = zeros(sampleN,gateCycleL+length(y));
smpl_pt = round(size(t_resp,2)/2); % take the mid point
s = zeros(sampleN,1);
if indep==0  % amp and phase not indep
    parfor iS=1:sampleN
        time = round(rand(1,N)*gateCycleL);
        if s_distr==0
            r = ones(1,N);
        else
            r = raylrnd(ones(1,N)/sqrt(2));
        end
        resp_temp = zeros(1,gateCycleL+length(y));
        for iN=1:N
            temp = [zeros(1,time(iN)),r(iN)*y,...
                    zeros(1,gateCycleL-time(iN))];
            resp_temp = resp_temp+temp;
        end
        resp_temp = abs(hilbert(resp_temp));
        %        resp(iS,:) = resp_temp;
        s(iS) = resp_temp(smpl_pt);
    end
else  % amp and phase indep
    parfor iS=1:sampleN
        time = round(rand(1,N)*gateCycleL);
        if s_distr==0
            r = ones(1,N);
        else
            r = raylrnd(ones(1,N)/sqrt(2));
        end
        resp_temp = zeros(1,gateCycleL+length(y));
        ph = cos(rand(1,N)*2*pi);  % random phase
        for iN=1:N
            temp = [zeros(1,time(iN)),r(iN)*y*sin(ph(iN)),...
                    zeros(1,gateCycleL-time(iN))];
            resp_temp = resp_temp+temp;
        end
        resp_temp = abs(hilbert(resp_temp));
        %        resp(iS,:) = resp_temp;
        s(iS) = resp_temp(smpl_pt);
    end    
end
save([result_dir,'/',result_file],...
     's','fs','sampleN','N','cycleNum','gateCycleNum','indep','s_distr');
