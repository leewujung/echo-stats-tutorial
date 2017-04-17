% 2014 09 01  batch calculate narrowband echopdf using Rayleigh
%             scatterers in order to compare with broadband
%             counterpart


%N = [10:10:90,100:100:3000]
N = [1800:100:3000];
ns = 1e5;
nb_freq = 50e3;
gate_len = 0.5;
sdir = '/mnt/storage/ECHO_STAT/20140901_nbechopdf_rayleigh_scatterers_1e5';
if ~exist(sdir,'dir')
    mkdir(sdir);
end

for iN=1:length(N)
    disp(['N=',num2str(N(iN))]);
    nbechopdf(N(iN),ns,gate_len,nb_freq,sdir,'rayleigh');
end

