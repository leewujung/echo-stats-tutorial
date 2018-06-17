% This code is for the creation of Fig. 10(b) in the tutorial.
% The equations used in the code below aer mathematically
% equivalent to eq.(34) and eq.(35) in the paper, although they are
% different in form. 
% Written by BAIK, Kyungmin. 10/18/206

ae=0.1; % radius of equal volume of sphere in m
Ep=0.1; % aspect ratio of spheroid

fes=ae/2; % backscattering amplitude for equal volume of sphere.

Fmax=fes*Ep^(-2/3);
Fmin=fes*Ep^(4/3);

N=1000;
Fran=linspace(1.01*Fmin,0.99*Fmax,N);

%PDF=((Ep^(2/3))*fes)./(2*(Fran.^(3/2))*sqrt(1-Ep^2).*sqrt(fes-(Ep^(2/3))*Fran)); % PDF for spheroid 3D
PDF=((Ep^(2/3))*fes)./(pi*Fran.*sqrt(fes-(Ep^(2/3))*Fran).*sqrt(Fran-(Ep^(4/3))*fes)); % 2D space

%PDF=PDF/trapz(Fran,PDF); % normalization to confirm integration of PDF over Fran is unity.

NT=1e4;

Ftotal=logspace(-5,2,NT);
PDFtotal=zeros(size(Ftotal));

for m=1:NT
    Rran=fliplr(Ftotal(m)./Fran); % Range of Rayleigh distribution
    PDFtotal(m)=trapz(Rran,exp(-0.5*(Rran.^2)).*fliplr(PDF));
end

% Sfac=trapz(Ftotal,PDFtotal);
% PDFtotal=PDFtotal/Sfac;

Nfac=sqrt(trapz(Ftotal,(Ftotal.^2).*PDFtotal)); % normalization factor for fs

loglog(Ftotal/Nfac,PDFtotal/trapz(Ftotal/Nfac,PDFtotal),'r','LineWidth',2)
axis([1e-3 1e2 1e-6 1e3])

hold on

Ep=0.2; % aspect ratio of spheroid

Fmax=fes*Ep^(-2/3);
Fmin=fes*Ep^(4/3);

Fran=linspace(1.01*Fmin,0.99*Fmax,N);

%PDF=((Ep^(2/3))*fes)./(2*(Fran.^(3/2))*sqrt(1-Ep^2).*sqrt(fes-(Ep^(2/3))*Fran)); % PDF for spheroid 3D
PDF=((Ep^(2/3))*fes)./(pi*Fran.*sqrt(fes-(Ep^(2/3))*Fran).*sqrt(Fran-(Ep^(4/3))*fes)); % 2D space

%PDF=PDF/trapz(Fran,PDF); % normalization to confirm integration of PDF over Fran is unity.

Ftotal=logspace(-5,2,NT);
PDFtotal=zeros(size(Ftotal));

for m=1:NT
    Rran=fliplr(Ftotal(m)./Fran); % Range of Rayleigh distribution
    PDFtotal(m)=trapz(Rran,exp(-0.5*(Rran.^2)).*fliplr(PDF));
end

Nfac=sqrt(trapz(Ftotal,(Ftotal.^2).*PDFtotal)); % normalization factor for fs
loglog(Ftotal/Nfac,PDFtotal/trapz(Ftotal/Nfac,PDFtotal),'g','LineWidth',2)

Ep=0.5; % aspect ratio of spheroid

Fmax=fes*Ep^(-2/3);
Fmin=fes*Ep^(4/3);

Fran=linspace(1.01*Fmin,0.99*Fmax,N);

%PDF=((Ep^(2/3))*fes)./(2*(Fran.^(3/2))*sqrt(1-Ep^2).*sqrt(fes-(Ep^(2/3))*Fran)); % PDF for spheroid 3D
PDF=((Ep^(2/3))*fes)./(pi*Fran.*sqrt(fes-(Ep^(2/3))*Fran).*sqrt(Fran-(Ep^(4/3))*fes)); % 2D space

%PDF=PDF/trapz(Fran,PDF); % normalization to confirm integration of PDF over Fran is unity.

Ftotal=logspace(-5,2,NT);
PDFtotal=zeros(size(Ftotal));

for m=1:NT
    Rran=fliplr(Ftotal(m)./Fran); % Range of Rayleigh distribution
    PDFtotal(m)=trapz(Rran,exp(-0.5*(Rran.^2)).*fliplr(PDF));
end

Nfac=sqrt(trapz(Ftotal,(Ftotal.^2).*PDFtotal)); % normalization factor for fs
loglog(Ftotal/Nfac,PDFtotal/trapz(Ftotal/Nfac,PDFtotal),'b','LineWidth',2)

Ftotal=logspace(-5,2,NT);
PDFtotal=Ftotal.*exp(-(Ftotal.^2)/2);

Nfac=sqrt(trapz(Ftotal,(Ftotal.^2).*PDFtotal)); % normalization factor for fs
loglog(Ftotal/Nfac,PDFtotal/trapz(Ftotal/Nfac,PDFtotal),'k','LineWidth',2)

legend('10:1','5:1','2:1','1:1 (sphere)')
axis([1e-3 1e2 1e-6 1e3])
%grid on
xlabel('$|{\it f_{ss.}}|/<|{\it f_{ss.}}|^{2}>^{1/2}$','interpreter','latex','FontSize',14)
ylabel('${\it p_{ss.}}(|{\it f_{ss.}}|)$','interpreter','latex','FontSize',14)