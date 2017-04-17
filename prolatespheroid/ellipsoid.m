function [pdf_x2, pdf_y2, cdf, pfa, theta, sigma_bs] = ellipsoid(e_ac)
% Calculates scattering amplitude of prolate spheroid according to angle.

% Setting constants
r = 0.1;
c = ((r^3)/(e_ac^2))^(1/3); 
e_ba = 1;
theta = linspace(0.0001, pi/2 - 0.0001, 1e6);

% Calculating scattering cross section
sigma_bs = (((c^2)*(e_ba^2))/4)*((sin(atan(e_ac./tan(theta)))./cos(theta)).^4);

% Converting to PDF
protopdf = sqrt(sigma_bs); % Converting from scattering cross section to scattering amplitude
pdfsize = 1e4;
pdf_y = zeros(1,pdfsize);
pdf_x = zeros(1,pdfsize);
inc2 = length(theta)/pdfsize;

parfor i = 1:pdfsize - 1
    y = protopdf(floor(inc2*i):floor(inc2*(i + 1)));
    x = theta(floor(inc2*i):floor(inc2*(i + 1)));
    if (length(y) > 1)
        pdf_y(i) = 1/(mean(abs(diff(y))));
        pdf_x(i) = mean(y);
    end
end
[pdf_x2 idx] = sort(pdf_x);
pdf_y2 = pdf_y(idx);

[pdf_x2, pdf_y2] = pdf_normalizer(pdf_x2(2:end), pdf_y2(2:end));

cdf = cumtrapz(pdf_x2, pdf_y2);
pfa = 1 - cdf;