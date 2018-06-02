% This code is for the creation of Fig. 14 in the tutorial.
% Written by BAIK, Kyungmin. 10/18/206

clear
t0=clock;

N=1e4;
b=logspace(-5,0,N);
theta=[0:.001:90]*pi/180;
PDF=zeros(size(b));
%ka=2*pi;
ka=44.251078712; 

% Beampattern PDF of a circular aperture 
for m=1:length(b)-1
    x=ka*sin(theta);
    y=4*(besselj(1,x)).^2-b(m)*(x.^2);
    thetaind=find(y(1:length(x)-1).*y(2:length(x))<0);
    thetaval=(theta(thetaind)+theta(thetaind+1))/2;
    PDF(m)=1/(2*pi*sqrt(b(m)))*sum(sin(thetaval)./(cos(thetaval).*abs(besselj(2,ka*sin(thetaval))))); % 2D
    %PDF(m)=1/(4*sqrt(b(m)))*sum(((sin(thetaval)).^2)./(cos(thetaval).*abs(besselj(2,ka*sin(thetaval))))); % 3D
end

figure(1)
Nfac=sqrt(trapz(b,(b.^2).*PDF)); % normalization factor for fs
bran=b/Nfac;
PDFnorm=PDF/trapz(bran,PDF);
loglog(bran,PDFnorm,'r','LineWidth',1)

hold on

figure(2)
% Beampattern PFA of a circular aperture
PFA=zeros(size(PDFnorm)-1);
for m=1:N-1
    PFA(m)=trapz(bran(m:N),PDFnorm(m:N));
end
loglog((bran(1:N-1)+bran(2:N))/2,PFA,'r','LineWidth',1)
hold on

ae=0.1; % radius of equal volume of sphere in m
fes=ae/2; % backscattering amplitude for equal volume of sphere.

% Combined with Rayleigh PDF
Nr=1e4;
Fray=logspace(-5,2,Nr);
Rran=Fray/fes;
PDFray=(Rran/fes).*exp(-0.5*(Rran.^2));

NT=1e4;
Fbeam=logspace(-5,0,NT);
%Ftotal=logspace(-20,log10(0.999999*Fmax),NT); % for smooth spheroid case
Ftotal=logspace(-9,2,NT);

% Combined with beampattern (Rayleigh scatterer)
PDFtotal=zeros(size(Ftotal));

for m=1:NT
    Bran=fliplr(Ftotal(m)./Fray); % Range of beampattern pdf
    %PDFtotal(m)=trapz(Bran,interp1(Ftotal,PDFbeam,Bran,'linear',0)./Bran.*fliplr(PDFray));
    PDFtotal(m)=trapz(Bran,interp1(Fbeam,PDF,Bran,'linear',0)./Bran.*fliplr(PDFray));
end

Nfac=sqrt(trapz(Ftotal,(Ftotal.^2).*PDFtotal)); % normalization factor for fs

% PDF of Rayleigh scatterer combined with beampattern PDF
figure (1)
loglog(Ftotal/Nfac,PDFtotal/trapz(Ftotal/Nfac,PDFtotal),'r','LineWidth',2)
axis([1e-2 1e2 1e-6 1e4])

% PFA of Rayleigh scatterer combined with beampattern PDF
figure(2)
PFA=zeros(size(PDFtotal)-1);
for m=1:NT-1
    PFA(m)=trapz(Ftotal(m:NT)/Nfac,PDFtotal(m:NT)/trapz(Ftotal/Nfac,PDFtotal));
end
loglog((Ftotal(1:NT-1)+Ftotal(2:NT))/(2*Nfac),PFA,'r','LineWidth',2)
axis([1e-2 1e2 1e-10 1e0])


PDFbeam=PDF;
Ep=0.1; % aspect ratio of spheroid

Fmax=fes*Ep^(-2/3);
Fmin=fes*Ep^(4/3);

Fran=logspace(log10(1.000001*Fmin),log10(0.999999*Fmax),N);

% PDF of spheroid in free-field
%PDF=((Ep^(2/3))*fes)./(2*(Fran.^(3/2))*sqrt(1-Ep^2).*sqrt(fes-(Ep^(2/3))*Fran)); % PDF for spheroid 3D
clear PDF
PDF=((Ep^(2/3))*fes)./(pi*Fran.*sqrt(fes-(Ep^(2/3))*Fran).*sqrt(Fran-(Ep^(4/3))*fes)); % 2D space

% This is for the calculation of smooth spheroid
PDFray=PDF; 
Fray=Fran;

Fbeam=logspace(-5,0,NT);
Ftotal=logspace(-9,2,NT);

% Combined with beampattern (smooth spheroid)
PDFtotal=zeros(size(Ftotal));
clear Bran
for m=1:NT
    Bran=fliplr(Ftotal(m)./Fray); % Range of beampattern pdf
    %PDFtotal(m)=trapz(Bran,interp1(Ftotal,PDFbeam,Bran,'linear',0)./Bran.*fliplr(PDFray));
    PDFtotal(m)=trapz(Bran,interp1(Fbeam,PDFbeam,Bran,'linear',0)./Bran.*fliplr(PDFray));
end

Nfac=sqrt(trapz(Ftotal,(Ftotal.^2).*PDFtotal)); % normalization factor for fs

% PDF of smooth spheroid combined with beampattern PDF 
figure(1)
loglog(Ftotal/Nfac,PDFtotal/trapz(Ftotal/Nfac,PDFtotal),'g','LineWidth',2)

% PFA of smooth spheroid combined with beampattern PDF
figure(2)
clear PFA
PFA=zeros(size(PDFtotal)-1);
for m=1:NT-1
    PFA(m)=trapz(Ftotal(m:NT)/Nfac,PDFtotal(m:NT)/trapz(Ftotal/Nfac,PDFtotal));
end
loglog((Ftotal(1:NT-1)+Ftotal(2:NT))/(2*Nfac),PFA,'g','LineWidth',2)


% Combined with Rayleigh PDF
clear Fray PDFray
Fray=logspace(-5,2,Nr);
PDFray=zeros(size(Fray));

for m=1:Nr
    Rran=fliplr(Fray(m)./Fran); % Range of Rayleigh distribution
    PDFray(m)=trapz(Rran,exp(-0.5*(Rran.^2)).*fliplr(PDF));
end

clear Fbeam Ftotal
Fbeam=logspace(-5,0,NT);
Ftotal=logspace(-9,2,NT);

% Combined with beampattern (Rough spheroid)
PDFtotal=zeros(size(Ftotal));
clear Bran
for m=1:NT
    Bran=fliplr(Ftotal(m)./Fray); % Range of beampattern pdf
    %PDFtotal(m)=trapz(Bran,interp1(Ftotal,PDFbeam,Bran,'linear',0)./Bran.*fliplr(PDFray));
    PDFtotal(m)=trapz(Bran,interp1(Fbeam,PDFbeam,Bran,'linear',0)./Bran.*fliplr(PDFray));
end

Nfac=sqrt(trapz(Ftotal,(Ftotal.^2).*PDFtotal)); % normalization factor for fs

% PDF of Rough spheroid combined with beampattern PDF in 2D
figure(1)
loglog(Ftotal/Nfac,PDFtotal/trapz(Ftotal/Nfac,PDFtotal),'r','LineWidth',2)

% PFA of Rough spheroid combined with beampattern PDF in 2D
figure(2)
clear PFA
PFA=zeros(size(PDFtotal)-1);
for m=1:NT-1
    PFA(m)=trapz(Ftotal(m:NT)/Nfac,PDFtotal(m:NT)/trapz(Ftotal/Nfac,PDFtotal));
end
loglog((Ftotal(1:NT-1)+Ftotal(2:NT))/(2*Nfac),PFA,'r','LineWidth',2)


eval(['disp(''total  ',num2str(etime(clock,t0)),' seconds elapsed'');'])
% legend('10:1','5:1','2:1','1:1')
% axis([1e-3 1e2 1e-6 1e3])
% grid on
% 
% xlabel('|{\it f_{ss}}|/<|{\it f_{ss}}|^2>^{1/2}','FontSize',14)
% ylabel('{\it p_{ss}}(|{\it f_{ss}}|)','FontSize',14)