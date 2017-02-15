% 2017 01 01  Rayleigh scatterer with and without beampattern, PDF & PFA

clear

%addpath '/mnt/storage/modeling_code_current/'
%addpath '/mnt/storage/analysis_code_current/'
addpath '~/Dropbox/0_CODE'/MATLAB/saveSameSize/

% base_path = '~/Desktop/echo_stat_figs';
base_path = '/Volumes/wjlee_apl_2/echo_stat_tutorial/echo_stat_figs/';

% Make save path
str = strsplit(mfilename('fullpath'),'/');
str = str{end};
save_path = fullfile(base_path,str);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Set param
X = load('fig_12_pb_ka_ka_num.mat');
ka = X.ka_3deg;
% ka = 2*pi;

pingnum_str = '1e7';
pingnum = eval(pingnum_str);

npt = 150;  % number of points for pe kde estimation

if ~exist(save_path,'dir')
    mkdir(save_path);
end

% N_all = 2:4;
% N_all = [1,2,5,10,50,100];
% N_all = [25,250,2500];
N_all = [1,10,100,1000];
v_rayl = 1/sqrt(2);

leg_str = cell(1,length(N_all)+1);
leg_str{1} = 'Rayleigh';
for iN=2:length(N_all)+1
    leg_str{iN} = sprintf('N=%d',N_all(iN-1));
end


if 0
    
for iN=1:length(N_all)
    Ns = N_all(iN);
    
    param.N = Ns;
    param.ka = ka;
    
    for iP = 1:pingnum
        phase = rand(1,Ns)*2*pi;
        amp = raylrnd(repmat(v_rayl,1,Ns));
        s = amp.*exp(1i*phase);
        
        % position in the beam
        count = 1;
        theta = zeros(1,Ns);
        while count <= Ns
            xx = rand(1);
            yy = rand(1);
            zz = rand(1);
            if sqrt(xx.^2+yy.^2+zz.^2)<1
                xy = sqrt(xx.^2+yy.^2);
                theta(count) = atan(xy./(zz));
                count = count +1;
            end
        end
        b = (2*besselj(1,ka*sin(theta))./(ka*sin(theta))).^2;
        
        % E=SB
        e = s.*b;
        env(iP) = abs(sum(e));
        
    end % pingnum
    
    file_save = sprintf('pnum_%s_ka%2.4f_N%04d.mat',...
        pingnum_str,ka,Ns);
    save([save_path,'/',file_save],'env','param'); 
end

end






% Plot: PDF WITH BEAMPATTERN
fig = figure;
xr = logspace(-3,log10(2000),500);  % standard
rayl = raylpdf(xr,1/sqrt(2));
loglog(xr,rayl,'k','linewidth',1);
hold on

for iN=1:length(N_all)
    
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(save_path,simu_file));
    %[x,p_x] = findEchoDist(E.env/sqrt(mean(E.env.^2)),npt);
    [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
    loglog(x,p_x,'-','linewidth',1);
    
end
xlabel('Normalized echo amplitude','fontsize',16);
ylabel('PDF','fontsize',16);
title(sprintf('ka=%2.4f, smplN=%s, with bp',...
    ka,pingnum_str),...
    'fontsize',18);
ll = legend(leg_str);
set(ll,'fontsize',18);
set(gca,'fontsize',14)
xlim([1e-3 1e2]);
ylim([1e-6 1e3]);

save_fname = sprintf('%s_ka%2.4f_smpl%s_pdf_bp1',...
    str,ka,pingnum_str);
saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
saveSameSize_100(fig,'file',[fullfile(save_path,save_fname),'.png'],...
    'format','png');


% Plot: PFA WITH BEAMPATTERN
fig = figure;
xr = logspace(-3,log10(2000),500);  % standard
rayl = raylpdf(xr,1/sqrt(2));
loglog(xr,rayl,'k','linewidth',1);
hold on

for iN=1:length(N_all)
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(save_path,simu_file));
    %[x,p_x] = findEchoDist(E.env/sqrt(mean(E.env.^2)),npt);
    [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
    cdf_x = cumtrapz(x,p_x);
    pfa_x = 1-cdf_x;
    loglog(x,pfa_x,'-','linewidth',1);
end
xlabel('Normalized echo amplitude','fontsize',16);
ylabel('PFA','fontsize',16);
title(sprintf('ka=%2.4f, smplN=%s, with bp',...
    ka,pingnum_str),...
    'fontsize',18);
ll = legend(leg_str);
set(ll,'fontsize',18);
set(gca,'fontsize',14)
xlim([1e-3 1e2]);
ylim([1e-6 1e3]);

save_fname = sprintf('%s_ka%2.4f_smpl%s_pfa_bp1',...
    str,ka,pingnum_str);
saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
saveSameSize_100(fig,'file',[fullfile(save_path,save_fname),'.png'],...
    'format','png');
