% Code to generate Figure 18 of the echo statistics tutorial.
% This figure shows PDF of magnitude of echo from 100 identical Rayleigh
% scatterers that are randomly and uniformly distributed in thin
% hemispherical shell 
%
% Author: Wu-Jung Lee | leewujung@gmail.com | APL-UW



clear
addpath './util_fcn'
base_path = './figs';

% Make save path
str = strsplit(mfilename('fullpath'),'/');
str = str{end};
save_path = fullfile(base_path,str);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Set param
pingnum_str = '1e7';
pingnum = eval(pingnum_str);
v_rayl = 1/sqrt(2);

npt = 80;  % number of points for pe kde estimation
N = 100;

X = load('./figs/figure_12/figure_12_ka_num.mat');
deg = [1,3,10,20];

% Set operation
mc_opt = 1;  % 0 - do not re-generate realizations
             % 1 - re-generate all realizations

% Monte Carlo simulation
if mc_opt
    for iD=1:length(deg)
        eval(['ka = X.ka_',num2str(deg(iD)),'deg;']);
        param.N = N;
        param.ka = ka;
        
        parfor iP = 1:pingnum
            phase = rand(1,N)*2*pi;
            amp = raylrnd(repmat(v_rayl,1,N));
            s = amp.*exp(1i*phase);
            
            % position in the beam
            u = unifrnd(0,1,1,sum(N));
            theta = acos(u);  % polar angle wrt beam axis
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
ll = legend('Rayleigh','1^o (0.0417)','3^o (0.375)',...
            '10^o (4.13)','20^o (16.0)');
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
saveSameSize(fig,'file',[fullfile(save_path,save_fname),'.png'],...
    'format','png');

