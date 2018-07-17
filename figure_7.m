% Code to generate Figure 7 of the echo statistics tutorial
%
% This code plots the Rice PDF on both a lin-lin and log-log scale.
%
% Author:  Kyungmin Baik | kbaik@kriss.re.kr | KRISS


clear all
close all

%Gamma=input('Enter the power signal-to-noise ratio: ');
X=logspace(-5,3,1e6);

Gamma=0;
PDF=2*(1+Gamma)*X.*exp(-((1+Gamma)*(X.^2)+Gamma)).*besseli(0,2*sqrt(Gamma*(1+Gamma))*X);
figure(1)
plot(X,PDF,'k-','LineWidth',2)
hold on
figure(2)
loglog(X,PDF,'k-','LineWidth',2)
hold on

Gamma=0.5;
PDF=2*(1+Gamma)*X.*exp(-((1+Gamma)*(X.^2)+Gamma)).*besseli(0,2*sqrt(Gamma*(1+Gamma))*X);
figure(1)
plot(X,PDF,'k-','LineWidth',1)
figure(2)
loglog(X,PDF,'k-','LineWidth',1)

Gamma=1;
PDF=2*(1+Gamma)*X.*exp(-((1+Gamma)*(X.^2)+Gamma)).*besseli(0,2*sqrt(Gamma*(1+Gamma))*X);
figure(1)
plot(X,PDF,'b-','LineWidth',1)
figure(2)
loglog(X,PDF,'b-','LineWidth',1)

Gamma=2;
PDF=2*(1+Gamma)*X.*exp(-((1+Gamma)*(X.^2)+Gamma)).*besseli(0,2*sqrt(Gamma*(1+Gamma))*X);
figure(1)
plot(X,PDF,'b-','LineWidth',2)
figure(2)
loglog(X,PDF,'b-','LineWidth',2)

Gamma=5;
PDF=2*(1+Gamma)*X.*exp(-((1+Gamma)*(X.^2)+Gamma)).*besseli(0,2*sqrt(Gamma*(1+Gamma))*X);
figure(1)
plot(X,PDF,'g-','LineWidth',2)
figure(2)
loglog(X,PDF,'g-','LineWidth',2)

Gamma=10;
PDF=2*(1+Gamma)*X.*exp(-((1+Gamma)*(X.^2)+Gamma)).*besseli(0,2*sqrt(Gamma*(1+Gamma))*X);
figure(1)
plot(X,PDF,'r-','LineWidth',2)
axis([0 3 0 2])
legend('\gamma=0 (Rayleigh)','\gamma=0.5','\gamma=1','\gamma=2','\gamma=5','\gamma=10')
figure(2)
loglog(X,PDF,'r-','LineWidth',2)
axis([1e-2 1e1 1e-6 1e2])
