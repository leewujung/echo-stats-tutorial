% Code to generate Figure 20 of the echo statistics tutorial.
% This figure shows PDF of echo magnitude from multiple Rayleigh scatterers
% in a split aggregation in which the larger scatterers are separated from
% the smaller scatterers as illustrated in Fig. 19a
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

A = 0.05;
M_all = [5,20];
Nw_all = 2500;
Ns_all = [25,250,2500];%,20,50,500];
v_rayl1 = 1/sqrt(2);
pingnum_str = '1e7';
pingnum = eval(pingnum_str);

npt = 120;  % number of points for pe kde estimation


% Set operation
mc_opt = 1;  % 0 - do not re-generate realizations
             % 1 - re-generate all realizations


% Monte Carlo simulation
if mc_opt
    for iM = 1%:length(M_all)
        v_rayl2 = M_all(iM)/sqrt(2);
        param.M = M_all(iM);
        disp(['M=',num2str(param.M)]);
        
        for iKA = 1:length(ka)
            disp(['ka=',num2str(ka(iKA))]);
            for iNw = 1:length(Nw_all)                
                for iNs = 1:length(Ns_all)
                    tic
                    ka_sl = ka(iKA);
                    Ns_sl = Ns_all(iNs);
                    Nw_sl = Nw_all(iNw);
                    param.ka = ka_sl;
                    param.Ns = Ns_sl;
                    param.Nw = Nw_all(iNw);
                    param.A  = A;

                    fprintf('Ns=%d, Nw=%d\n',Ns_sl,Nw_sl);
                    
                    pingnum1 = pingnum*(1-param.A);
                    pingnum2 = pingnum*param.A;
                    env1 = zeros(1,length(pingnum1));
                    env2 = zeros(1,length(pingnum2));
                    
                    parfor iP = 1:pingnum1
                        % SCATTERER 1
                        % before beampattern
                        v_rayl = v_rayl1;
                        phase = rand(1,Nw_sl)*2*pi;
                        amp = raylrnd(repmat(v_rayl,1,Nw_sl));
                        s1 = amp.*exp(1i*phase);
                        
                        % position in the beam
                        u = unifrnd(0,1,1,sum(Nw_sl));
                        theta = acos(u);  % polar angle wrt beam axis
                        b1 = (2*besselj(1,ka_sl*sin(theta))./(ka_sl*sin(theta))).^2;
                        
                        % E=SB
                        e1 = s1.*b1;
                        env1(iP) = abs(sum(e1));
                    end
                    
                    parfor iP = 1:pingnum2
                        % SCATTERER 2
                        % before beampattern
                        v_rayl = v_rayl2;
                        phase = rand(1,Ns_sl)*2*pi;
                        amp = raylrnd(repmat(v_rayl,1,Ns_sl));
                        s2 = amp.*exp(1i*phase);
                        
                        % position in the beam
                        u = unifrnd(0,1,1,sum(Ns_sl));
                        theta = acos(u);  % polar angle wrt beam axis
                        b2 = (2*besselj(1,ka_sl*sin(theta))./(ka_sl*sin(theta))).^2;
                        
                        % E=SB
                        e2 = s2.*b2;
                        env2(iP) = abs(sum(e2));
                        
                    end
                    env = [env1,env2];
                    file_save = sprintf('pnum_%s_ka%2.4f_A%2.2f_M%02d_Nw%04d_Ns%04d.mat',...
                        pingnum_str,ka(iKA),A,param.M, ...
                        param.Nw,param.Ns);
                    save([save_path,'/',file_save],'env','param');
                    
                    toc
                end  % Ns
            end  % Nw
        end  % ka
    end  % scale
end % if re-run simulation




% Plot and cmp
for iM=1:length(M_all)
    for iNw=1:length(Nw_all)
        for iNs=1:length(Ns_all)
            fprintf('Ns=%d, Nw=%d\n',Ns_all(iNs),Nw_all(iNw));
            save_fname = sprintf('%s_A%2.2f_M%02d_Nw%04d_Ns%04d_smpl%s',...
                str,A,M_all(iM),Nw_all(iNw),Ns_all(iNs),pingnum_str);
            
            fig = figure;
            xr = logspace(-3,log10(2000),500);  % standard
            rayl = raylpdf(xr,1/sqrt(2));
            hr = loglog(xr,rayl,'k','linewidth',2);
            hold on
            
            simu_file = sprintf('pnum_%s_ka%2.4f_A%2.2f_M%02d_Nw%04d_Ns%04d.mat',...
                pingnum_str,ka,A,M_all(iM),Nw_all(iNw),Ns_all(iNs));
            E = load(fullfile(save_path,simu_file));
            [p_x,x] = findEchoDist_kde(E.env/sqrt(mean(E.env.^2)),npt);
            hh = loglog(x,p_x,'r','linewidth',2);
            
            %     title(sprintf('A=%2.2f ,Ns=%d, Nw=%d, smplN=%s',...
            %                   A,Ns_all(iNs),Nw_all(iNw),pingnum_str),...
            %           'fontsize',18);
            set(gca,'fontsize',16)
            xlabel('$\tilde{e}/<\tilde{e}^2>^{1/2}$','Interpreter','LaTex','fontsize',24);
            ylabel('$p_e(\tilde{e}/<\tilde{e}^2>^{1/2})$','Interpreter','LaTex','fontsize',24);
            switch iNs
                case 1
                    ll = legend(hh,'Ns = 25 (0.0937)');
                case 2
                    ll = legend(hh,'Ns = 250 (0.937)');
                case 3
                    ll = legend(hh,'Ns = 2500 (9.37)');
            end
            set(ll,'fontsize',22);
            xlim([1e-3 1e2]);
            ylim([1e-6 1e3]);
            
            saveas(fig,[fullfile(save_path,save_fname),'.fig'],'fig');
            saveSameSize(fig,'file',[fullfile(save_path,save_fname),'.png'],...
                'format','png');
            
        end
    end
end
