% 2017 01 01  Echo pdf of 100 Rayleigh scatterers observed using different
%             beamwidth
% 2017 04 12  Update figure legend, axis labels, and curve style


addpath '~/code_matlab_dn/saveSameSize/'

% base_path = '~/Desktop/echo_stat_figs';
% base_path = '/Volumes/wjlee_apl_2/echo_stat_tutorial/echo_stat_figs/';
base_path = '/home/wu-jung/internal_2tb/echo_stat_tutorial/echo_stat_figs';

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
loglog(xr,rayl,'k','linewidth',2);
hold on

for iD=1:length(deg)
    eval(['ka = X.ka_',num2str(deg(iD)),'deg;']);
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d.mat',...
        pingnum_str,ka,N);
    E = load(fullfile(save_path,simu_file));
    %[x,p_x] = findEchoDist(E.env/sqrt(mean(E.env.^2)),npt);
    [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
    switch iD
        case 1
            loglog(x,p_x,'r-','linewidth',2);   
        case 2
            loglog(x,p_x,'g-','linewidth',2);   
        case 3
            loglog(x,p_x,'b-','linewidth',2);   
        case 4
            loglog(x,p_x,'b-','linewidth',1);   
    end
    clear E
end
% title(sprintf('N=%d, smplN=%s',...
%     N,pingnum_str),...
%     'fontsize',18);
ll = legend('Rayleigh','1^o','3^o','10^o','20^o');
set(ll,'fontsize',18,'location','southwest');
new_pos = ll.Position;
new_pos(1) = 0.3;
set(ll,'Position',new_pos);
text(3e-3,3e2,'100 Rayleigh scatterers','fontsize',24);
set(gca,'fontsize',16)
xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
ylabel('$p_e(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
xlim([1e-3 1e2]);
ylim([1e-6 1e3]);

save_fname = sprintf('%s_smpl%s',...
    str,pingnum_str);
saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
saveSameSize_100(fig,'file',[fullfile(save_path,save_fname),'.png'],...
    'format','png');

