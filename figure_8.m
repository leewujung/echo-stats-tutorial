% This code plots the PDF of K-distribution in normalized quantity
% Written by BAIK, Kyungmin. 11/14/2016

clear

Alpha=input('Enter the value of alpha: '); % Shape parameter
X=logspace(-5,3,1e6);

PDF=4*sqrt(Alpha)/gamma(Alpha)*((sqrt(Alpha)*X).^Alpha).*besselk(Alpha-1,2*sqrt(Alpha)*X);

%figure
%plot(X,PDF,'k-','LineWidth',2)
%axis([0 3 0 2])
%grid on
%xlabel('|A|/<|A|^2>^{1/2} ','FontSize',12)
%ylabel('Probability Density Function','FontSize',12)
%title('K-Distribution','FontSize',12)

% 
% figure
loglog(X,PDF,'k-','LineWidth',2)
%axis([1e-2 1e2 1e-6 1e4])
axis([1e-2 2e2 1e-4 1e2])
grid on
xlabel('|A|/<|A|^2>^{1/2} ','FontSize',12)
ylabel('Probability Density Function','FontSize',12)
% title('K-Distribution','FontSize',12)

%xlabel('|e_{sys}|/<|e_{sys}|^2>^{1/2} ','FontSize',12)
%ylabel('Probability Density Function','FontSize',12)

%xlabel('A_{norm} ','FontSize',12)
%ylabel('p_{R,norm}(A_{norm})','FontSize',12)