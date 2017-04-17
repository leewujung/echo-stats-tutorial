function pm = gather_bbechopdf_curves(N,sdir,ns,glen,fpre,npt,pmx)
%
% Compile broadband echo pdf models into one matrix for
% curve-fitting
%
% INPUT 
%   N      number of scatterers in the gate
%   sdir   directory of model files
%   ns     number of samples
%   glen   lenght of gate [m]
%   fpre   prefix of model filename
%   npt    number of points for each model curve using KDE
%   pmx    x-axis for the model curves
% OUTPUT
%   pm     model curves obtained by KDE
%
% Wu-Jung Lee  2014/04/27  revive from old code


tail = sprintf('_sampleN%d_gateLen%2.1f_freqDepBP.mat',...
               ns,glen);
%pmx = logspace(-5,100,1000);  % pre-set x-axis of model curves

pm = zeros(length(pmx),length(N));
for iN=1:length(N)
    fname = [fpre,'_N_',num2str(N(iN)),'_r_1',tail];
    A = load([sdir,'/',fname]);
    [pm_tmp,pmx_tmp] = findEchoDist_kde(A.s/sqrt(mean(A.s.^2)),npt);
    pm_tmp(isnan(pm_tmp)==1) = [];
    pmx_tmp(isnan(pm_tmp)==1) = [];    
    pm(:,iN) = interp1(pmx_tmp,pm_tmp,pmx);
end


