% This code creates 6 common PDFs in one figure in linear or log scale

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
loglog(X,PDF,'r-','LineWidth',2)
axis([1e-2 5e2 1e-3 1e3])

subplot(322)
Alpha=0.5;
PDF=4*sqrt(Alpha)/gamma(Alpha)*((sqrt(Alpha)*X).^Alpha).*besselk(Alpha-1,2*sqrt(Alpha)*X);
loglog(X,PDF,'k-','LineWidth',2)
axis([1e-2 5e2 1e-3 1e3])

subplot(323)
Gamma=0.5;
PDF=Gamma*((2/Gamma*gamma(2/Gamma))^(Gamma/2))*(X.^(Gamma-1)).*exp(-((2/Gamma*gamma(2/Gamma))^(Gamma/2))*(X.^Gamma));
loglog(X,PDF,'c-','LineWidth',1)
axis([1e-2 5e2 1e-3 1e3])

subplot(324)
Sig=0.2;
PDF=1./(Sig*sqrt(2*pi)*X).*exp(-((log(X)+Sig^2).^2)/(2*(Sig^2)));
loglog(X,PDF,'m-','LineWidth',1)
axis([1e-2 5e2 1e-3 1e3])

subplot(325)
m=0.5;
PDF=2*(m^m)/gamma(m)*(X.^(2*m-1)).*exp(-m*(X.^2));
loglog(X,PDF,'m-','LineWidth',1)
axis([1e-2 5e2 1e-3 1e3])

subplot(326)
Gamma=-2;
if Gamma==0
    PDF=sqrt(2/((1-2*Gamma)*(1-Gamma)))*exp(-X*sqrt(2/((1-2*Gamma)*(1-Gamma))));
else
    PDF=sqrt(2/((1-2*Gamma)*(1-Gamma)))*((1+Gamma*X*sqrt(2/((1-2*Gamma)*(1-Gamma)))).^(-(1+Gamma)/Gamma));
end
loglog(X,PDF,'m-','LineWidth',1)
axis([1e-2 5e2 1e-3 1e3])

