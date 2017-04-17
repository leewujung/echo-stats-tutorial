% 2013 07 30  Generate truncated normal distribution samples

function sout = normrnd_truncated(mu,std,nstd,ns,s)

stmp = normrnd(mu,std,1,ns);

if any(stmp>mu+nstd*std | stmp<mu-nstd*std)
    stmp(stmp>(mu+nstd*std) | stmp<(mu-nstd*std)) = [];
    stmp = normrnd_truncated(mu,std,nstd,ns-length(stmp),stmp);
end

sout = [s,stmp];