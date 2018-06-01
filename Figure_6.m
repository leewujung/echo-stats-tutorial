% This code plots the PDF of magnitude of signal due to sums of N simple sinusoids with identical
% magnitude and random phases.
% Written by BAIK, Kyungmin. 10/25/2016

clear

N=100; % number of identical random variables
Ns=1e7; % number of samples
Ampv=zeros(Ns,1);
Noisamp=0.3;
%Noisamp=0;
%%Ampr=zeros(Ns,1);
%%Ampi=zeros(Ns,1);
%%Std=1; % standard deviation of amplitude

for m=1:N
    Pha=2*pi*rand(Ns,1); % random phase within 0 and 2*pi
    Nois=-Noisamp+2*Noisamp*rand(Ns,1);
    %%Phar=2*pi*rand(Ns,1); % random phase within 0 and 2*pi
    %%Phai=2*pi*rand(Ns,1); % random phase within 0 and 2*pi
    %Famp=1-Std+2*Std*rand(Ns,1); % random amplitude between 1-Std and 1+Std
    Famp=1; % random amplitude between 1-Std and 1+Std
    %Ampv=Ampv+Famp*exp(1i*Pha); % amplitude of each phase component
    Ampv=Ampv+Nois+Famp*exp(1i*Pha); % amplitude of each phase component, noise is added
    %%Ampr=Ampr+Famp*cos(Phar); % amplitude of each phase component
    %%Ampi=Ampi+Famp*sin(Phai); % amplitude of each phase component
    
    
end

%Amp=sqrt(Ampr.^2+Ampi.^2);
Amp=abs(Ampv);
%Amp=Ampv;

Nbin=3e2; % number of magnitude bin
%Nbin=120; % number of magnitude bin



%Edges=linspace(0,Famp*N,Nbin+1); % edges for bins of amplitude in log-scale

Edges=linspace(0,Famp*N+Noisamp,Nbin+1); % edges for bins of amplitude in log-scale

%%Edges=logspace(log10(min(Amp)),log10(max(Amp)),Nbin+1); % edges for bins of amplitude in log-scale

Df=diff(Edges); % Delta Amp

PDFfs=zeros(Nbin,2);

[AA,BB]=hist(Amp,Edges);

%AA=histc(Amp,Edges);
PDFfs(:,1)=sqrt(BB(1:Nbin).*BB(2:Nbin+1));
%%PDFfs(:,1)=BB'/(sqrt(sum(BB.^2)/Nbin));
PDFfs(:,2)=(AA(1:Nbin)./Df)';

Sfac=trapz(PDFfs(:,1),PDFfs(:,2));
PDFfs(:,2)=PDFfs(:,2)/Sfac;

Nfac=sqrt(trapz(PDFfs(:,1),(PDFfs(:,1).^2).*PDFfs(:,2))); % normalization factor for fs

%%loglog(PDFfs(:,1)/Nfac,PDFfs(:,2)/trapz(PDFfs(:,1)/Nfac,PDFfs(:,2)),'k.','LineWidth',1)
plot(PDFfs(:,1)/Nfac,PDFfs(:,2)/trapz(PDFfs(:,1)/Nfac,PDFfs(:,2)),'k-','LineWidth',2)
%%plot(PDFfs(1:40:Nbin,1)/Nfac,PDFfs(1:40:Nbin,2)/trapz(PDFfs(1:40:Nbin,1)/Nfac,PDFfs(1:40:Nbin,2)),'bo','MarkerSize',6)
%loglog(PDFfs(:,1),PDFfs(:,2),'r.','LineWidth',1)
%axis([1e-3 1e2 1e-6 1e3])
grid on