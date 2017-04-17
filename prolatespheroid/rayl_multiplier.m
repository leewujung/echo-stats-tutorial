function [x_pdf, y_pdf, cdf, pfa] = rayl_multiplier(x, fx)
% Written by Saurav Bhatia, WHOI, spring 2011.
% Finds PDF of the product v of random variables x and y. Assumes y is
% Rayleigh distributed.
% Integrates over x. 

v = logspace(-6, 4, 3e4);
fv = zeros(1, length(v));

xn = x(x ~= 0);
fxn = fx(x ~= 0);

parfor ii = 1:length(v) % Evaluates fv at various values of v.
    fyr = raylpdf(v(ii)./xn, 1/sqrt(2)); % Determines fy(v/x) for use in the integral.
    acc = fxn.*fyr./(xn); % Evaluates the term within the integral at all values of x.
    % acc(xn < v(ii)) = zeros(1, length(acc(xn < v(ii)))); % Setting limits of integration
    fv(ii) = trapz(xn, acc); % Evaluates the integral at one value of v.
end

[x_pdf y_pdf] = pdf_normalizer(v, fv);

cdf = cumtrapz(x_pdf, y_pdf);
pfa = 1 - cdf;
