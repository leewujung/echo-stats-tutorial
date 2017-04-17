function [pdf_x pdf_y theta bp] = cyl_bp_num(ka)
% Numerically calculates PDF of the beampattern of cylindrical piston transducer given k*a value

theta = logspace(-5, log10(pi/2 - 1e-8), 2e6);

% Calculating beampattern for piston transducer
nump = 2*besselj(1, ka*sin(theta));
denp = 2*pi*sin(theta);
bp = (nump./denp).^2;

% Numerically calculating PDF from beampattern
pdfsize = 2e4;
inc = max(abs(bp))/pdfsize;
temp_y = zeros(1,pdfsize);
temp_x = zeros(1,pdfsize);

parfor i = 1:pdfsize
    idx = (bp >= (i - 1)*inc) & (bp < i*inc);
    y = bp(idx);
    if (length(y) > 1)
        num = 2;
        temp_y(i) = num/(pi*mean(abs(diff(y))));
        temp_x(i) = mean(y);
    end
end

[temp_x2 ind] = sort(temp_x);
temp_y2 = temp_y(ind);
[pdf_x pdf_y] = pdf_normalizer(temp_x2, temp_y2);

% Debug plots
% close all
% figure(1)
% semilogy(theta, bp, '.')
% title(sprintf('Beampattern intensity, ka = %upi', ka/pi))
% grid on
% set(gca, 'XMinorGrid', 'Off', 'YMinorGrid','Off');
% axis([1e-4 pi/2 1e-10 10])
% figure(2)
% loglog(pdf_x, pdf_y, '.')
% title(sprintf('Beampattern PDF, ka = %upi', ka/pi))
% grid on
% set(gca, 'XMinorGrid', 'Off', 'YMinorGrid','Off');