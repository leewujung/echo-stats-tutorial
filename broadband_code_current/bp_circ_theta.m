function bp = bp_circ_theta(theta,freq,a)
% Plot beampattern as a function of thera and freq

k = 2*pi*freq'/1500;
bp = (2*besselj(1,k*a*sin(theta))./(k*a*sin(theta))).^2;
figure;
imagesc(theta/pi*180,freq/1e3,10*log10(abs(bp)));
xlabel('\theta {degree})','fontsize',14);
ylabel('Frequency (kHz)','fontsize',14);
set(gca,'fontsize',12);
colorbar
