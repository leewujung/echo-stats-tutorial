% 2013 07 29  Generate samples from discrete distribution
%             according to Handbook of Monte Carlo Methods (Kroese,
%             Taimre, and Botev, 2011), info obtained from Matlab central
%             file exchange

function [s,idx] = discrete_rnd(x,p,ns)

if size(p,2)~=1
    p = p';
end

[~,idx] = histc(rand(1,ns),[0;cumsum(p)]); 
s = x(idx);
if size(s,2)~=1
    s=s';
end