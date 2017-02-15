% 2017 01 01  Echo pdf of 100 Rayleigh scatterers observed using different
%             beamwidth

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

pingnum_str = '1e7';
pingnum = eval(pingnum_str);
v_rayl = 1/sqrt(2);

npt = 80;  % number of points for pe kde estimation
N = 100;

if ~exist(save_path,'dir')
    mkdir(save_path);
end

X = load('fig_12_pb_ka_ka_num.mat');
deg = [1,3,10,20];


if 0
    
for iD=1:length(deg)
    eval(['ka = X.ka_',num2str(deg(iD)),'deg;']);
    param.N = N;
    param.ka = ka;
    
    for iP = 1:pingnum
        phase = rand(1,N)*2*pi;
        amp = raylrnd(repmat(v_rayl,1,N));
        s = amp.*exp(1i*phase);
        
        % position in the beam
        count = 1;
        theta = zeros(1,N);
        while count <= N
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
        pingnum_str,ka,N);
    save([save_path,'/',file_save],'env','param');
end

end




% Plot: PDF
fig = figure;
xr = logspace(-3,log10(2000),500);  % standard
rayl = raylpdf(xr,1/sqrt(2));
loglog(xr,rayl,'k','linewidth',1);
hold on

for iD=1:length(deg)
    eval(['ka = X.ka_',num2str(deg(iD)),'deg;']);
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d.mat',...
        pingnum_str,ka,N);
    E = load(fullfile(save_path,simu_file));
    %[x,p_x] = findEchoDist(E.env/sqrt(mean(E.env.^2)),npt);
    [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
    loglog(x,p_x,'-','linewidth',1);
    
end
xlabel('Normalized echo amplitude','fontsize',16);
ylabel('PDF','fontsize',16);
title(sprintf('N=%d, smplN=%s',...
    N,pingnum_str),...
    'fontsize',18);
ll = legend('Rayleigh','1^o','3^o','10^o','20^o');
set(ll,'fontsize',18);
set(gca,'fontsize',14)
xlim([1e-3 1e2]);
ylim([1e-6 1e3]);

save_fname = sprintf('%s_smpl%s',...
    str,pingnum_str);
saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
saveSameSize_100(fig,'file',[fullfile(save_path,save_fname),'.png'],...
    'format','png');

