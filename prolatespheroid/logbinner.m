function [x_pdf y_pdf] = logbinner(data, n, x_min)
% Divides data into n logarithmically spaced bins with an offset of
% 10^(x_min)

edges = logspace(log10(min(data)) + x_min, log10(max(data)), n);
y_pdf = histc(data, edges);
y_pdf = y_pdf(1:end-1); % Removing last bin, since it only counts values that match edges(end).

bins = sqrt(edges(1:end-1).*edges(2:end));
[x_pdf, y_pdf] = pdf_normalizer(bins, (y_pdf)./diff(edges));