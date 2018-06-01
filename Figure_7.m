% This code plots the Rice PDF in normalized quantity
% Written by BAIK, Kyungmin. 11/14/2016

clear

Gamma=input('Enter the power signal-to-noise ratio: ');
X=logspace(-5,3,1e6);

PDF=2*(1+Gamma)*X.*exp(-((1+Gamma)*(X.^2)+Gamma)).*besseli(0,2*sqrt(Gamma*(1+Gamma))*X);

%figure
plot(X,PDF,'r-','LineWidth',2)
% dashline(X,PDF,'k-',2,1,2,1,'LineWidth',1)
% axis([0 3 0 2])
% grid on
% xlabel('|A|/<|A|^2>^{1/2} ','FontSize',12)
% ylabel('Probability Density Function','FontSize',12)
% title('Rice Distribution','FontSize',12)


% figure
%loglog(X,PDF,'r-','LineWidth',2)
%dashline(X,PDF,2,1,2,1,'k-','LineWidth',1)
%axis([1e-2 1e2 1e-6 1e4])
%grid on
%xlabel('|A|/<|A|^2>^{1/2} ','FontSize',12)
%ylabel('Probability Density Function','FontSize',12)
%title('Rice Distribution','FontSize',12)



%xlabel('|e_{sys}|/<|e_{sys}|^2>^{1/2} ','FontSize',12)
%ylabel('Probability Density Function','FontSize',12)

%xlabel('A_{norm} ','FontSize',12)
%ylabel('p_{C,norm}(A_{norm})','FontSize',12)