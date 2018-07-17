% Code to generate Figure 10a of the echo statistics tutorial
%
% This code plots the PDF of the magnitude of the backscattering amplitude
% for the randomly oriented smooth prolate spheroid. The equations used in
% the code are mathematically equivalent to Eqs. (33) and (34) in the paper,
% although they are different in form.
%
% Author:  Kyungmin Baik | kbaik@kriss.re.kr | KRISS


clear
ae=0.1; % radius of equal volume of sphere in m
Ep=0.1; % aspect ratio of spheroid

fes=ae/2; % backscattering amplitude for equal volume of sphere.

Fmax=fes*Ep^(-2/3);
Fmin=fes*Ep^(4/3);

N=10000;
%Fran=linspace(1.001*Fmin,0.999*Fmax,N);
Fran=linspace(1.0000001*Fmin,0.9999999*Fmax,N);

%PDF=((Ep^(2/3))*fes)./(2*(Fran.^(3/2))*sqrt(1-Ep^2).*sqrt(fes-(Ep^(2/3))*Fran)); % 3D space
PDF=((Ep^(2/3))*fes)./(pi*Fran.*sqrt(fes-(Ep^(2/3))*Fran).*sqrt(Fran-(Ep^(4/3))*fes)); % 2D space

% Sfac=trapz(Fran,PDF);
% PDF=PDF/Sfac;

Nfac=sqrt(trapz(Fran,(Fran.^2).*PDF)); % normalization factor for fs

loglog(Fran/Nfac,PDF/trapz(Fran/Nfac,PDF),'r','LineWidth',2)


hold on

Ep=0.2; % aspect ratio of spheroid

Fmax=fes*Ep^(-2/3);
Fmin=fes*Ep^(4/3);

Fran=linspace(1.0000001*Fmin,0.9999999*Fmax,N);

%PDF=((Ep^(2/3))*fes)./(2*(Fran.^(3/2))*sqrt(1-Ep^2).*sqrt(fes-(Ep^(2/3))*Fran)); % 3D space
PDF=((Ep^(2/3))*fes)./(pi*Fran.*sqrt(fes-(Ep^(2/3))*Fran).*sqrt(Fran-(Ep^(4/3))*fes)); % 2D space

Nfac=sqrt(trapz(Fran,(Fran.^2).*PDF)); % normalization factor for fs

loglog(Fran/Nfac,PDF/trapz(Fran/Nfac,PDF),'g','LineWidth',2)

Ep=0.5; % aspect ratio of spheroid

Fmax=fes*Ep^(-2/3);
Fmin=fes*Ep^(4/3);

Fran=linspace(1.0000001*Fmin,0.9999999*Fmax,N);

%PDF=((Ep^(2/3))*fes)./(2*(Fran.^(3/2))*sqrt(1-Ep^2).*sqrt(fes-(Ep^(2/3))*Fran)); % 3D space
PDF=((Ep^(2/3))*fes)./(pi*Fran.*sqrt(fes-(Ep^(2/3))*Fran).*sqrt(Fran-(Ep^(4/3))*fes)); % 2D space

Nfac=sqrt(trapz(Fran,(Fran.^2).*PDF)); % normalization factor for fs

loglog(Fran/Nfac,PDF/trapz(Fran/Nfac,PDF),'b','LineWidth',2)

%legend('10:1','5:1','2:1')
%grid on

xlabel('$|{\it f_{ss}}|/<|{\it f_{ss}}|^2>^{1/2}$','interpreter','latex','FontSize',14)
ylabel('${\it p_{ss}}(|{\it f_{ss}}|)$','interpreter','latex','FontSize',14)
