% 2012 05 15  Compile all modeling results into one matrix to facilitate
%             curve-fitting
% 2012 09 25  Update this using kde density estimates
% 2012 10 17  Option to add noise

function [pm,pmx,N,s_all] = get_model_curves_kde(BB_MODEL_DIR,N,npt,fpost,varargin)
%npt = 200;
% varargin{1}  threshold
% varargin{2}  ratio for Rayleigh noise
%folder = '/mnt/storage/ECHO_STAT/20120814_bb_1e4smpl_2';
%folder = '/mnt/storage/ECHO_STAT/20121007_bb_5e4smpl';
%folder = '/mnt/storage/ECHO_STAT/20121026_bb_smpl1e4_noise';
folder = BB_MODEL_DIR;

sampleN = 5e4;
%gateLdist = 0.5;
%N=[10,20,50,100:200:900,2000];
%N=[10:10:100,150:50:1000];
addpath '/mnt/storage/analysis_code_current';

fname = 'tx_sys';
p = zeros(length(N),npt);
pm_temp = cell(1,length(N));
pmx_temp = cell(1,length(N));
s_all = cell(1,length(N));
%s_all = zeros(length(N),sampleN);
for iN=1:length(N)
    disp(['N=',num2str(N(iN))]);
    %fload = [folder,'/',fname,'_','N',num2str(N(iN)),...
    %         '_sampleN',num2str(sampleN),...
    %         '_gateLdist',num2str(gateLdist),'_freqDepBP.mat'];
    %fload = [folder,'/',fname,'_','N',num2str(N(iN)),...
    %         '_sampleN',num2str(sampleN),...
    %         '_gateLen',num2str(gateLdist),'_freqDepBP.mat'];
    fload = [folder,'/',fname,'_','N',num2str(N(iN)),...
             '_sampleN',num2str(sampleN),'_',fpost];
    if ~isempty(varargin)
        if length(varargin)==2
            tmp = load(fload,'s','resp_x');
            tmplen = length(tmp.s);
            if tmplen>1e4
                seg = 5;
                seglen = floor(tmplen/seg);
                if mod(tmplen,seg)~=0;
                    seg = seg+1;
                end
                s = zeros(1,tmplen);
                tic
                for iS=1:seg
                    disp(['seg=',num2str(iS)]);
                    sidx = (1:seglen)+(iS-1)*seglen;
                    ss = zeros(1,tmplen);
                    ss(sidx) = 1;
                    ss_logical = logical(ss);
                    sidx = find((ss_logical & logical(ones(1,tmplen)))==1);
                    resp_x = tmp.resp_x(sidx,:);
                    n = normrnd(0,1/sqrt(2)/varargin{2},size(resp_x,1),size(resp_x,2));
                    resp_x = resp_x+n;
                    resp_x_env = abs(hilbert(resp_x'))';
                    mid_idx = floor((size(resp_x,2)+1)/2);
                    smpl = resp_x_env(:,mid_idx);
                    %s = [s,resp_x_env(:,mid_idx)];
                    s(sidx) = smpl;
                end
                toc
            end

            % note: varargin{2} = frms/nrms
            %n = normrnd(0,1/sqrt(2)/varargin{2},size(tmp.resp_x,1),size(tmp.resp_x,2));
            %resp_x = tmp.resp_x+n;
            %resp_x_env = abs(hilbert(resp_x'))';
            %mid_idx = floor((size(resp_x,2)+1)/2);
            %s = resp_x_env(:,mid_idx);
            if ~isempty(varargin{1})
                s(s<varargin{1}) = [];
            end
        elseif length(varargin)==1
            tmp = load(fload,'s');
            s = tmp.s;
            if ~isempty(varargin{1})
                s(s<varargin{1}) = [];
            end
        end
    else
        tmp = load(fload,'s');
        s = tmp.s;
    end

    %    [x,ptemp] = findEchoDist(s,'default',0,npt);
    s = s/sqrt(mean(s.^2));  % normalization

    [pm_temp{iN},pmx_temp{iN},~] = findEchoDist_kde(s,npt);
    %    nanidx = find(isnan(pm_temp{iN})==1);  % take out NaN points
    %    pm_temp{iN}(nanidx) = [];
    %    pmx_temp{iN}(nanidx) = [];
    pmx_range(iN) = pmx_temp{iN}(end)-pmx_temp{iN}(1);
    %    s_all(iN,:) = s;
    s_all{iN} = s;
end

[~,xidx] = max(pmx_range);
pmx = pmx_temp{xidx};  % the output x-axis
pm = nan(length(pmx),length(N));
for iN=1:length(N)
    pm(:,iN) = interp1(pmx_temp{iN},pm_temp{iN},pmx);
end

% add Rayleigh distribution
N(end+1) = Inf;
pm(:,length(N)) = raylpdf(pmx,1/sqrt(2));

%{
figure;
loglog(pmx,pm(:,end),'k');
hold on
for iN=1:length(N)-1
    loglog(pmx_temp{iN},pm_temp{iN});
end
%}