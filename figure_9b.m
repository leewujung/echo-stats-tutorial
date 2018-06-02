% This code is for the creation of Fig. 9(b) in the tutorial.
% Written by BAIK, Kyungmin. 10/18/206

clear
ae=0.1; % radius of equal volume of sphere in m
Ep=0.1; % aspect ratio of spheroid

fes=ae/2; % backscattering amplitude for equal volume of sphere.
Gamma=linspace(0,90,1000)*pi/180;
Amps=0.5*ae*(Ep^(4/3))./((Ep*cos(Gamma)).^2+(sin(Gamma)).^2); % amplitude of spheroid
semilogy(Gamma*180/pi,Amps,'r','LineWidth',2)

hold on
Ep=0.2; % aspect ratio of spheroid
Amps=0.5*ae*(Ep^(4/3))./((Ep*cos(Gamma)).^2+(sin(Gamma)).^2); % amplitude of spheroid
semilogy(Gamma*180/pi,Amps,'g','LineWidth',2)

Ep=0.5; % aspect ratio of spheroid
Amps=0.5*ae*(Ep^(4/3))./((Ep*cos(Gamma)).^2+(sin(Gamma)).^2); % amplitude of spheroid
semilogy(Gamma*180/pi,Amps,'b','LineWidth',2)

Ep=1; % aspect ratio of spheroid
Amps=0.5*ae*(Ep^(4/3))./((Ep*cos(Gamma)).^2+(sin(Gamma)).^2); % amplitude of spheroid
semilogy(Gamma*180/pi,Amps,'k','LineWidth',2)

xlabel('\beta','FontSize',14)
ylabel('|{\it f_{ss}}|','FontSize',14)
grid on

axis([0 90 1e-3 1e0])
legend('10:1','5:1','2:1','1:1 (sphere)')