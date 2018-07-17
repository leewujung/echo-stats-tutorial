% Code to generate Figure A1 of the echo statistics tutorial
%
% This code plots 6 common (“generic”) PDFs in one figure on a lin-lin scale.
% The Rayleigh PDF is given in each panel for comparison.
%
% Author:  Kyungmin Baik | kbaik@kriss.re.kr | KRISS


clear

FigHandle = figure(1);
set(FigHandle, 'Position', [100, 100, 850, 900]);
subplot('Position', [3/48, 11/16, 10/24, 0.275])
subplot('Position', [27/48, 11/16, 10/24, 0.275])
subplot('Position', [3/48, 5.75/16, 10/24, 0.275])
subplot('Position', [27/48, 5.75/16, 10/24, 0.275])
subplot('Position', [3/48, 0.5/16, 10/24, 0.275])
subplot('Position', [27/48, 0.5/16, 10/24, 0.275])

X=logspace(-5,3,1e6);

subplot(321)
Gamma=0;
PDF=2*(1+Gamma)*X.*exp(-((1+Gamma)*(X.^2)+Gamma)).*besseli(0,2*sqrt(Gamma*(1+Gamma))*X);
plot(X,PDF,'k-','LineWidth',2)
hold on

Gamma=0.5;
PDF=2*(1+Gamma)*X.*exp(-((1+Gamma)*(X.^2)+Gamma)).*besseli(0,2*sqrt(Gamma*(1+Gamma))*X);
plot(X,PDF,'k-','LineWidth',1)

Gamma=1;
PDF=2*(1+Gamma)*X.*exp(-((1+Gamma)*(X.^2)+Gamma)).*besseli(0,2*sqrt(Gamma*(1+Gamma))*X);
plot(X,PDF,'b-','LineWidth',2)

Gamma=2;
PDF=2*(1+Gamma)*X.*exp(-((1+Gamma)*(X.^2)+Gamma)).*besseli(0,2*sqrt(Gamma*(1+Gamma))*X);
plot(X,PDF,'b-','LineWidth',2)

Gamma=5;
PDF=2*(1+Gamma)*X.*exp(-((1+Gamma)*(X.^2)+Gamma)).*besseli(0,2*sqrt(Gamma*(1+Gamma))*X);
plot(X,PDF,'g-','LineWidth',2)

Gamma=10;
PDF=2*(1+Gamma)*X.*exp(-((1+Gamma)*(X.^2)+Gamma)).*besseli(0,2*sqrt(Gamma*(1+Gamma))*X);
plot(X,PDF,'r-','LineWidth',2)
axis([0 3 0 2])
legend('\gamma=0 (Rayleigh)','\gamma=0.5','\gamma=1','\gamma=2','\gamma=5','\gamma=10')

subplot(322)
Alpha=0.5;
PDF=4*sqrt(Alpha)/gamma(Alpha)*((sqrt(Alpha)*X).^Alpha).*besselk(Alpha-1,2*sqrt(Alpha)*X);
plot(X,PDF,'r--','LineWidth',2)
hold on

Alpha=1;
PDF=4*sqrt(Alpha)/gamma(Alpha)*((sqrt(Alpha)*X).^Alpha).*besselk(Alpha-1,2*sqrt(Alpha)*X);
plot(X,PDF,'g--','LineWidth',2)

Alpha=2;
PDF=4*sqrt(Alpha)/gamma(Alpha)*((sqrt(Alpha)*X).^Alpha).*besselk(Alpha-1,2*sqrt(Alpha)*X);
plot(X,PDF,'b--','LineWidth',2)

Alpha=5;
PDF=4*sqrt(Alpha)/gamma(Alpha)*((sqrt(Alpha)*X).^Alpha).*besselk(Alpha-1,2*sqrt(Alpha)*X);
plot(X,PDF,'b--','LineWidth',1)

Alpha=100;
PDF=4*sqrt(Alpha)/gamma(Alpha)*((sqrt(Alpha)*X).^Alpha).*besselk(Alpha-1,2*sqrt(Alpha)*X);
plot(X,PDF,'b','LineWidth',2)
legend('\alpha_K=0.5 (exponential)','\alpha_K=1','\alpha_K=2','\alpha_K=5','\alpha_K=\infty (Rayleigh)')
axis([0 3 0 2])

subplot(323)
nu=0.5;
PDF=nu*((2/nu*gamma(2/nu))^(nu/2))*(X.^(nu-1)).*exp(-((2/nu*gamma(2/nu))^(nu/2))*(X.^nu));
plot(X,PDF,'r--','LineWidth',2)
hold on

nu=1;
PDF=nu*((2/nu*gamma(2/nu))^(nu/2))*(X.^(nu-1)).*exp(-((2/nu*gamma(2/nu))^(nu/2))*(X.^nu));
plot(X,PDF,'g--','LineWidth',2)

nu=1.5;
PDF=nu*((2/nu*gamma(2/nu))^(nu/2))*(X.^(nu-1)).*exp(-((2/nu*gamma(2/nu))^(nu/2))*(X.^nu));
plot(X,PDF,'b--','LineWidth',2)

nu=2;
PDF=nu*((2/nu*gamma(2/nu))^(nu/2))*(X.^(nu-1)).*exp(-((2/nu*gamma(2/nu))^(nu/2))*(X.^nu));
plot(X,PDF,'k','LineWidth',2)

nu=3;
PDF=nu*((2/nu*gamma(2/nu))^(nu/2))*(X.^(nu-1)).*exp(-((2/nu*gamma(2/nu))^(nu/2))*(X.^nu));
plot(X,PDF,'g','LineWidth',2)

nu=5;
PDF=nu*((2/nu*gamma(2/nu))^(nu/2))*(X.^(nu-1)).*exp(-((2/nu*gamma(2/nu))^(nu/2))*(X.^nu));
plot(X,PDF,'r','LineWidth',2)
legend('\nu=0.5','\nu=1 (exponential)','\nu=1.5','\nu=2 (Rayleigh)','\nu=3','\nu=5')
axis([0 3 0 2.5])

subplot(324)
Sig=0.2;
PDF=1./(Sig*sqrt(2*pi)*X).*exp(-((log(X)+Sig^2).^2)/(2*(Sig^2)));
plot(X,PDF,'r','LineWidth',2)
hold on

Sig=0.5;
PDF=1./(Sig*sqrt(2*pi)*X).*exp(-((log(X)+Sig^2).^2)/(2*(Sig^2)));
plot(X,PDF,'g-','LineWidth',2)

Sig=1;
PDF=1./(Sig*sqrt(2*pi)*X).*exp(-((log(X)+Sig^2).^2)/(2*(Sig^2)));
plot(X,PDF,'b--','LineWidth',2)

Sig=1.5;
PDF=1./(Sig*sqrt(2*pi)*X).*exp(-((log(X)+Sig^2).^2)/(2*(Sig^2)));
plot(X,PDF,'g--','LineWidth',2)

Sig=3;
PDF=1./(Sig*sqrt(2*pi)*X).*exp(-((log(X)+Sig^2).^2)/(2*(Sig^2)));
plot(X,PDF,'r--','LineWidth',2)

% Add Rayleigh
plot(X,2*X.*exp(-X.^2),'k','LineWidth',2)
legend('\sigma_{LN}=0.2','\sigma_{LN}=0.5','\sigma_{LN}=1','\sigma_{LN}=1.5','\sigma_{LN}=3','Rayleigh')
axis([0 3 0 2.5])

subplot(325)
m=0.5;
PDF=2*(m^m)/gamma(m)*(X.^(2*m-1)).*exp(-m*(X.^2));
plot(X,PDF,'r--','LineWidth',2)
hold on

m=1;
PDF=2*(m^m)/gamma(m)*(X.^(2*m-1)).*exp(-m*(X.^2));
plot(X,PDF,'k','LineWidth',2)

m=2;
PDF=2*(m^m)/gamma(m)*(X.^(2*m-1)).*exp(-m*(X.^2));
plot(X,PDF,'b','LineWidth',2)

m=5;
PDF=2*(m^m)/gamma(m)*(X.^(2*m-1)).*exp(-m*(X.^2));
plot(X,PDF,'g','LineWidth',2)

m=10;
PDF=2*(m^m)/gamma(m)*(X.^(2*m-1)).*exp(-m*(X.^2));
plot(X,PDF,'r','LineWidth',2)
legend('m=0.5','m=1 (Rayleigh)','m=2','m=5','m=10')
axis([0 3 0 3])

subplot(326)
Rho=-2;
PDF=sqrt(2/((1-2*Rho)*(1-Rho)))*((1+Rho*X*sqrt(2/((1-2*Rho)*(1-Rho)))).^(-(1+Rho)/Rho));
plot(X(PDF>0),PDF(PDF>0),'r','LineWidth',2)
hold on

Rho=-1;
PDF=sqrt(2/((1-2*Rho)*(1-Rho)))*((1+Rho*X*sqrt(2/((1-2*Rho)*(1-Rho)))).^(-(1+Rho)/Rho));
plot(X(PDF>0),PDF(PDF>0),'g','LineWidth',2)

Rho=-0.5;
PDF=sqrt(2/((1-2*Rho)*(1-Rho)))*((1+Rho*X*sqrt(2/((1-2*Rho)*(1-Rho)))).^(-(1+Rho)/Rho));
plot(X(PDF>0),PDF(PDF>0),'b--','LineWidth',2)

Rho=0;
PDF=sqrt(2/((1-2*Rho)*(1-Rho)))*exp(-X*sqrt(2/((1-2*Rho)*(1-Rho))));
plot(X(PDF>0),PDF(PDF>0),'g--','LineWidth',2)

Rho=1/3;
PDF=sqrt(2/((1-2*Rho)*(1-Rho)))*((1+Rho*X*sqrt(2/((1-2*Rho)*(1-Rho)))).^(-(1+Rho)/Rho));
plot(X(PDF>0),PDF(PDF>0),'r--','LineWidth',2)

% Add Rayleigh
plot(X,2*X.*exp(-X.^2),'k','LineWidth',2)
legend('\rho=-2','\rho=-1','\rho=-1/2','\rho=0','\rho=1/3','Rayleigh')
axis([0 3 0 3])
