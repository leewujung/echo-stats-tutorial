% Code to generate Figure 8 of the echo statistics tutorial
%
% This code plots the K PDF on both a lin-lin and log-log scale.
% The equation used in the code is mathematically equivalent
% to Eq. (29) in the paper, although it is slightly different in
% form.
%
% Author:  Kyungmin Baik | kbaik@kriss.re.kr | KRISS


clear all
close all

X=logspace(-5,3,1e6);
%Alpha=input('Enter the value of alpha: '); % Shape parameter
Alpha=0.5;
PDF=4*sqrt(Alpha)/gamma(Alpha)*((sqrt(Alpha)*X).^Alpha).*besselk(Alpha-1,2*sqrt(Alpha)*X);
figure(1)
plot(X,PDF,'r--','LineWidth',2)
hold on
figure(2)
loglog(X,PDF,'r--','LineWidth',2)
hold on

Alpha=1;
PDF=4*sqrt(Alpha)/gamma(Alpha)*((sqrt(Alpha)*X).^Alpha).*besselk(Alpha-1,2*sqrt(Alpha)*X);
figure(1)
plot(X,PDF,'g--','LineWidth',2)
figure(2)
loglog(X,PDF,'g--','LineWidth',2)

Alpha=2;
PDF=4*sqrt(Alpha)/gamma(Alpha)*((sqrt(Alpha)*X).^Alpha).*besselk(Alpha-1,2*sqrt(Alpha)*X);
figure(1)
plot(X,PDF,'b--','LineWidth',2)
figure(2)
loglog(X,PDF,'b--','LineWidth',2)

Alpha=5;
PDF=4*sqrt(Alpha)/gamma(Alpha)*((sqrt(Alpha)*X).^Alpha).*besselk(Alpha-1,2*sqrt(Alpha)*X);
figure(1)
plot(X,PDF,'b--','LineWidth',1)
figure(2)
loglog(X,PDF,'b--','LineWidth',1)

Alpha=100; % sufficiently large value for the case of alpha_k=infinite
PDF=4*sqrt(Alpha)/gamma(Alpha)*((sqrt(Alpha)*X).^Alpha).*besselk(Alpha-1,2*sqrt(Alpha)*X);
figure(1)
plot(X,PDF,'k-','LineWidth',2)
axis([0 3 0 2])
legend('\alpha_K=0.5 (exponential)','\alpha_K=1','\alpha_K=2','\alpha_K=5','\alpha_K=\infty (Rayleigh)')
figure(2)
loglog(X,PDF,'k-','LineWidth',2)
axis([1e-2 1e1 1e-4 1e1])
