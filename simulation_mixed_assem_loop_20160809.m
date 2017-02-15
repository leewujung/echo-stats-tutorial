% 2010 12 12  numerical simulation of Nscat scatterers in the beam
% 2016 08 09  simulation for echo stat tutorial


%addpath '/mnt/storage/broadband_code_current/'
%addpath '/mnt/storage/analysis_code_current/'

save_path_pre = '~/wjl/echo_stat_figs';

% Make save path
str = strsplit(mfilename('fullpath'),'/');
str = str{end};


setnum = 2;
pingnum_str = '1e5';
pingnum = eval(pingnum_str);
Nw = 100;
Ns = [1,10,50];
v_rayl1 = 1/sqrt(2);
%scale_all = [5:5:60];
scale_all = 20;

for iSet=1:setnum
    disp(['Set=',num2str(iSet)]);

    save_path = fullfile(save_path_pre,sprintf('%s_set%d',str,iSet));
    if ~exist(save_path,'dir')
        mkdir(save_path);
    end

    for iS = 1:length(scale_all)
        disp(['str=',num2str(scale_all(iS))]);
        v_rayl2 = scale_all(iS)/sqrt(2);
        param.M = scale_all(iS);

        ka = [5,8,14];

        for iKA = 1:length(ka)
            disp(['ka=',num2str(ka(iKA)),'pi']);
            
            for iN = 1:length(Ns)
                tic
                env = zeros(1,length(pingnum));
                ka_sl = ka(iKA);
                Ns_sl = Ns(iN);
                param.ka = ka_sl*pi;
                param.Ns = Ns_sl;
                param.Nw = Nw;

                disp(['Ns=',num2str(Ns_sl)]);
                parfor iP = 1:pingnum
                    % SCATTERER 1
                    % before beampattern
                    v_rayl = v_rayl1;
                    phase = rand(1,Nw)*2*pi;
                    amp = raylrnd(repmat(v_rayl,1,Nw));
                    s1 = amp.*exp(1i*phase);
                    
                    % position in the beam
                    count = 1;
                    theta = zeros(1,Nw);
                    while count <= Nw
                        xx = rand(1);
                        yy = rand(1);
                        zz = rand(1);
                        if sqrt(xx.^2+yy.^2+zz.^2)<1
                            xy = sqrt(xx.^2+yy.^2);
                            theta(count) = atan(xy./(zz));
                            count = count +1;
                        end
                    end
                    b1 = (2*besselj(1,ka_sl*pi*sin(theta))./(ka_sl*pi*sin(theta))).^2;
                    
                    % E=SB
                    e1 = s1.*b1;
                    
                    % SCATTERER 2
                    % before beampattern
                    v_rayl = v_rayl2;
                    phase = rand(1,Ns_sl)*2*pi;
                    amp = raylrnd(repmat(v_rayl,1,Ns_sl));
                    s2 = amp.*exp(1i*phase);
                    
                    % position in the beam
                    count = 1;
                    theta = zeros(1,Ns_sl);
                    while count <= Ns_sl
                        xx = rand(1);
                        yy = rand(1);
                        zz = rand(1);
                        if sqrt(xx.^2+yy.^2+zz.^2)<1
                            xy = sqrt(xx.^2+yy.^2);
                            theta(count) = atan(xy./(zz));
                            count = count +1;
                        end
                    end
                    b2 = (2*besselj(1,ka_sl*pi*sin(theta))./(ka_sl*pi*sin(theta))).^2;
                    
                    % E=SB
                    e2 = s2.*b2;
                    
                    % SUMMATION
                    env(iP) = abs(sum([e1,e2]));
                    
                end % pingnum

                file_save = sprintf('pnum_%s_ka%dpi_M%d_Nw%d_Ns%d.mat',...
                                    pingnum_str,ka(iKA),param.M, ...
                                    param.Nw,param.Ns);

                save([save_path,'/',file_save],'env','param');

                toc
            end % Ns

        end % ka

    end % scale


end % iSet