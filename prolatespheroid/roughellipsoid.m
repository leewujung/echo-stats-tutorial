% close all
function [pdf_x3, pdf_y3, cdf, pfa] = roughellipsoid(e_ac)
% Calculates scattering amplitude of prolate spheroid according to angle.

[pdf_x2, pdf_y2] = ellipsoid(e_ac);

[pdf_x3, pdf_y3, cdf, pfa] = rayl_multiplier(pdf_x2, pdf_y2);
