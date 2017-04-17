function [xout, yout] = pdf_normalizer(x, y)
% Normalizes PDFs in X- and Y- directions
rmsval = sqrt((trapz(x, y.*x.*x))/trapz(x, y));

xout = x/rmsval;

yout = y/trapz(xout, y);