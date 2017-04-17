function [x_pdf, y_pdf, cdf, pfa] = pdf_multiplier3(xb, yb, x2, y2)
% Written by Saurav Bhatia, WHOI, summer 2009.
% Combines  PDFs assuming xb, yb is beampattern PDF and y2 is the
% target response across all values of x2. Uses equation (15) 
% from paper by Chu, Stanton.

x2 = x2(x2 ~= 0);
y2 = y2(x2 ~= 0);

v = logspace(-5, 2, 5e2); % The output PDF is calculated on this set of points
y_pdf = zeros(1, length(v));

parfor ii = 1:length(v)-1
    y2_temp = interp1(x2, y2, v(ii)./xb, 'spline', 0);
    acc = yb.*y2_temp./xb;
    y_pdf(ii) = trapz(xb, acc);
end

[x_pdf, y_pdf] = pdf_normalizer(v, y_pdf);

cdf = cumtrapz(x_pdf, y_pdf);
pfa = 1 - cdf;
