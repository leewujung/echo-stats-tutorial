% 2012 12 02  Generate echo pdf from rough prolate spheroid, with noise

function [pm,pmx,pm0,pmx0] =  prosph_2D_noisy_pdf(ar,fnrms,smpln)

addpath /mnt/storage/prolatespheroid/

%ar = 5;

for iAR = 1:length(ar)
    [dmx0(:,iAR), dm0(:,iAR), s0] =...
        prosph_2D_simulation(1/ar(iAR),smpln,1);
    [pm0(:,iAR),pmx0(:,iAR),~] = findEchoDist_kde(s0/sqrt(mean(s0.^2)),1e2);
    
    [dmx(:,iAR), dm(:,iAR), s] =...
        prosph_2D_simulation(1/ar(iAR),smpln,1,fnrms);
    [pm(:,iAR),pmx(:,iAR),~] = findEchoDist_kde(s/sqrt(mean(s.^2)),1e2);
end

pmx_new = logspace(log10(min(min(pmx))),log10(max(max(pmx))),300);
pmx0_new = logspace(log10(min(min(pmx0))),log10(max(max(pmx0))),300);
for iAR = 1:length(ar)
    pm_new(:,iAR) = interp1(pmx(:,iAR),pm(:,iAR),pmx_new);
    pm0_new(:,iAR) = interp1(pmx0(:,iAR),pm0(:,iAR),pmx0_new);
end
pmx = pmx_new;
pmx0 = pmx0_new;
pm = pm_new;
pm0 = pm0_new;