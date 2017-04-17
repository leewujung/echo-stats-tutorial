% 2012 02 24  find zero-crossing points for the xcorr results

%tmp = diff(sign(x(smpl_pt)/max(x)));
tmp = diff(sign(data));
idxp = find(tmp>0.5);  % - to +
idxn = find(tmp<-0.5); % + to -
idx = sort([idxp,idxn]);