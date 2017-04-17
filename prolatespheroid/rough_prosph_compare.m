close all
load('data/tsdata.mat')
load('data/pdf_prolatespheroid.mat')
c = 1; 
e_ba = 1;
e_ac = 0.1;
theta = linspace(0.00001, pi/2 - 0.0001, 500000);
beta = theta(theta<pi/4);
[x_beampattern, y_beampattern, x_in, y_in] = cyl_bp_num(2*pi);

[x_10_1, y_10_1, theta_10_1, sigma_10_1] = roughellipsoid(0.1);
[x_5_1, y_5_1, theta_5_1, sigma_5_1] = roughellipsoid(0.2);
[x_3_1, y_3_1, theta_3_1, sigma_3_1] = roughellipsoid(1/3);
[x_2_1, y_2_1, theta_2_1, sigma_2_1] = roughellipsoid(0.5);

x_r = logspace(-5, 2, 10000);
y_r = raylpdf(x_r, 0.7072);

figure
loglog(x_10_1, y_10_1, 'k', 'LineWidth', 1.2)
grid on
hold on
loglog(x_5_1, y_5_1, 'k-.', 'LineWidth', 1.2)
% loglog(x_3_1, y_3_1, 'k-.')
loglog(x_2_1, y_2_1, 'k--', 'LineWidth', 1.2)
loglog(x_r, y_r, 'k:', 'LineWidth', 1.2)
axis([1e-5 1e2 1e-5 1e2])
legend('10:1', '5:1', '2:1', '1:1 (Rayleigh');
set(gca, 'XMinorGrid', 'Off', 'YMinorGrid','Off');
ylabel('Probability Density Function')
xlabel('Normalized Echo Amplitude')
title('Rough prolate spheroid echo PDF, no beampattern effects')

fprintf('Calculating prolate spheroid echo 10:1\n');
tic;
[x_10_1_bp, y_10_1_bp, cdf_10_1, pfa_10_1] = pdf_multiplier2(x_beampattern, y_beampattern, x_10_1, y_10_1);
time(1) = toc;
fprintf('Calculations complete in %g seconds\n', time(1));

fprintf('Calculating prolate spheroid echo 5:1\n');
tic;
[x_5_1_bp, y_5_1_bp, cdf_5_1, pfa_5_1] = pdf_multiplier2(x_beampattern, y_beampattern, x_5_1, y_5_1);
time(2) = toc;
fprintf('Calculations complete in %g seconds\n', time(2));

fprintf('Calculating prolate spheroid echo 3:1\n');
tic;
[x_3_1_bp, y_3_1_bp, cdf_3_1, pfa_3_1] = pdf_multiplier2(x_beampattern, y_beampattern, x_3_1, y_3_1);
time(3) = toc;
fprintf('Calculations complete in %g seconds\n', time(3));

fprintf('Calculating prolate spheroid echo 2:1\n');
tic;
[x_2_1_bp, y_2_1_bp, cdf_2_1, pfa_2_1] = pdf_multiplier2(x_beampattern, y_beampattern, x_2_1, y_2_1);
time(4) = toc;
fprintf('Calculations complete in %g seconds\n', time(4));

fprintf('Calculating Rayleigh echo\n');
tic;
[xr_bp, yr_bp, cdf_r, pfa_r] = pdf_multiplier2(x_beampattern, y_beampattern, x_r, y_r);
time(5) = toc;
fprintf('Calculations complete in %g seconds\n', time(5));

figure 
loglog(x_10_1_bp, y_10_1_bp, 'k')
grid on
hold on
loglog(x_5_1_bp, y_5_1_bp, 'k--')
loglog(x_3_1_bp, y_3_1_bp, 'k:')
loglog(x_2_1_bp, y_2_1_bp, 'k-.')
loglog(xr_bp, yr_bp, 'k', 'LineWidth', 2)
axis([1e-3 10 1e-7 100])
legend('10:1', '5:1', '3:1', '2:1', '1:1');
set(gca, 'XMinorGrid', 'Off', 'YMinorGrid','Off');
ylabel('Probability Density Function')
xlabel('Normalized Echo Amplitude')
title('Rough prolate spheroid echo PDF with beampattern effects')

figure 
loglog(x_10_1_bp, pfa_10_1)
grid on
hold on
loglog(x_5_1_bp, pfa_5_1, 'k')
loglog(x_3_1_bp, pfa_3_1, 'Color', [0.6 0 0])
loglog(x_2_1_bp, pfa_2_1, 'r')
loglog(xr_bp, pfa_r, 'g')
axis([1e-2 10 1e-7 10])
legend('10:1', '5:1', '3:1', '2:1', 'Rayleigh');
set(gca, 'XMinorGrid', 'Off', 'YMinorGrid','Off');
ylabel('Probability of False Alarm')
xlabel('Normalized Echo Amplitude')
title('Rough prolate spheroid echo PFA with beampattern effects')

% matlabpool close
