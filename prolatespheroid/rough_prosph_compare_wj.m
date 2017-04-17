%close all
%load('data/tsdata.mat')
%load('data/pdf_prolatespheroid.mat')
c = 1; 
e_ba = 1;
e_ac = 0.1;
theta = linspace(0.00001, pi/2 - 0.0001, 500000);
beta = theta(theta<pi/4);
[x_beampattern, y_beampattern, x_in, y_in] = cyl_bp_num(2*pi);

ar = 1:0.1:10;
for iAR=2:length(ar)
    [x, y(iAR,:),~,~] = roughellipsoid(1/ar(iAR));
end

y(1,:) = raylpdf(x,1/sqrt(2));
