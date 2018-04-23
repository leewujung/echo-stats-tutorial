% Code to generate Figure 15 of the echo statistics tutorial
% This figure shows Distributions of magnitude of echo in backscatter
% direction from N point scatterers randomly and uniformly distributed in
% a thin hemispherical shell  
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
X = load('./figs/figure_12/figure_12_ka_num.mat');
ka = X.ka_3deg;

pingnum_str = '1e8';
pingnum = eval(pingnum_str);

N_all = [1,10,100,1000];
v_rayl = 1/sqrt(2);

% Set operation
mc_opt = 1;  % 0 - do not re-generate realizations
             % 1 - re-generate all realizations

% Monte Carlo simulation
if mc_opt
    for iN=1:length(N_all)
        Ns = N_all(iN);
        
        param.N = Ns;
        param.ka = ka;
        
        for iP = 1:pingnum
            phase = rand(1,Ns)*2*pi;
            amp = ones(1,Ns);
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
            b_bp0 = 1;
            b_bp1 = (2*besselj(1,ka*sin(theta))./(ka*sin(theta))).^2;
            
            % E=SB
            e_bp0 = s.*b_bp0;
            e_bp1 = s.*b_bp1;
            env_bp0(iP) = abs(sum(e_bp0));
            env_bp1(iP) = abs(sum(e_bp1));
            
        end % pingnum
        
        env = env_bp0;
        file_save = sprintf('pnum_%s_ka%2.4f_N%04d_bp0.mat',...
            pingnum_str,ka,Ns);
        save([save_path,'/',file_save],'env','param');
        
        env = env_bp1;
        file_save = sprintf('pnum_%s_ka%2.4f_N%04d_bp1.mat',...
            pingnum_str,ka,Ns);
        save([save_path,'/',file_save],'env','param');
        
    end
end



% Plot: PDF WITH BEAMPATTERN
N_all = [1,10,100,1000];
fig = figure;
xr = logspace(-3,log10(2000),500);  % standard
rayl = raylpdf(xr,1/sqrt(2));
loglog(xr,rayl,'k','linewidth',2);
hold on
for iN=1:length(N_all)    
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_bp1.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(save_path,simu_file));
    [x,p_x] = findEchoDist(E.env/sqrt(mean(E.env.^2)),1200);
    switch iN
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
% title(sprintf('ka=%2.4f, smplN=%s, with bp',...
%     ka,pingnum_str),...
%     'fontsize',18);
ll = legend('Rayleigh','N=1 (0.00375)','N=10 (0.0375)',...
    'N=100 (0.375)','N=1000 (3.75)',...
    'location','southwest');
set(ll,'fontsize',18);
set(gca,'fontsize',16)
xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
ylabel('$p_e(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
xlim([1e-3 1e2]);
ylim([1e-6 1e3]);

save_fname = sprintf('%s_smpl%s_ka%2.4f_pdf_bp1',...
    str,pingnum_str,ka);
saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
saveSameSize_100(fig,'file',[fullfile(save_path,save_fname),'.png'],...
    'format','png');


% Plot: PFA WITH BEAMPATTERN
N_all = [1,10,100,1000];
fig = figure;
xr = logspace(-3,log10(2000),500);  % standard
rayl = raylpdf(xr,1/sqrt(2));
cdf_rayl = cumtrapz(xr,rayl);
pfa_rayl = 1-cdf_rayl;
loglog(xr,pfa_rayl,'k','linewidth',2);
hold on
for iN=1:length(N_all)    
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_bp1.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(save_path,simu_file));
    [x,p_x] = findEchoDist(E.env/sqrt(mean(E.env.^2)),600);
    cdf_x = cumtrapz(x,p_x);
    pfa_x = 1-cdf_x;
    switch iN
        case 1
            loglog(x,pfa_x,'r-','linewidth',2);   
        case 2
            loglog(x,pfa_x,'g-','linewidth',2);   
        case 3
            loglog(x,pfa_x,'b-','linewidth',2);   
        case 4
            loglog(x,pfa_x,'b-','linewidth',1);   
    end
    clear E
end
% title(sprintf('ka=%2.4f, smplN=%s, with bp',...
%     ka,pingnum_str),...
%     'fontsize',18);
ll = legend('Rayleigh','N=1 (0.00375)','N=10 (0.0375)',...
    'N=100 (0.375)','N=1000 (3.75)',...
    'location','southwest');
set(ll,'fontsize',18);
set(gca,'fontsize',16)
xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
ylabel('$PFA(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
xlim([1e-3 1e2]);
ylim([1e-4 1e1]);

save_fname = sprintf('%s_smpl%s_ka%2.4f_pfa_bp1',...
    str,pingnum_str,ka);
saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
saveSameSize_100(fig,'file',[fullfile(save_path,save_fname),'.png'],...
    'format','png');


% Plot: PDF NO BEAMPATTERN
N_all = 1:4;

fig = figure;
xr = logspace(-3,log10(2000),500);  % standard
rayl = raylpdf(xr,1/sqrt(2));
loglog(xr,rayl,'k','linewidth',2);
hold on
plot([1,1],[10e-8,1],'r-','linewidth',2);  % N=1 curve
for iN=2:length(N_all)    
    simu_file = sprintf('pnum_%s_ka%2.4f_N%04d_bp0.mat',...
        pingnum_str,ka,N_all(iN));
    E = load(fullfile(save_path,simu_file));
    [x,p_x] = findEchoDist(E.env/sqrt(mean(E.env.^2)),500);
    switch iN
        case 2
            loglog(x,p_x,'g-','linewidth',2);   
        case 3
            loglog(x,p_x,'b-','linewidth',2);   
        case 4
            loglog(x,p_x,'b-','linewidth',1);   
    end
    clear E
end
% title(sprintf('ka=%2.4f, smplN=%s, no bp',...
%     ka,pingnum_str),...
%     'fontsize',18);
ll = legend('Rayleigh','N=1','N=2','N=3','N=4',...
    'location','southwest');
set(ll,'fontsize',18);
new_pos = ll.Position;
new_pos(1) = 0.25;
set(ll,'Position',new_pos);
set(gca,'fontsize',16);
set(gca,'xscale','log','yscale','log');  % log scale
xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
ylabel('$p_e(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
xlim([1e-3 1e2]);
ylim([1e-6 1e3]);

save_fname = sprintf('%s_smpl%s_ka%2.4f_pdf_bp0',...
    str,pingnum_str,ka);
saveas(fig,[fullfile(save_path,save_fname),'_log.fig'],'fig');
saveSameSize_100(fig,'file',[fullfile(save_path,save_fname),'_log.png'],...
    'format','png');

