% 2014 04 27  Figure out where the code was left in August 2013
% 2014 05 26  Run more fish scatterer models for comparison
% 2014 09 04  Run bbechopdf using prolate spheroid scatterers with
%             realistic aspect ratio

mix_r = 1;
ns = 5e4;
gate_len = 0.5;
save_opt = 1;
tx = 3;
sdir = '/mnt/storage/ECHO_STAT/20140904_bbechopdf_prolate_sphroid_realistic_ar_length_variation';
if ~exist(sdir,'dir')
    mkdir(sdir);
end


N = [10:10:90,100:100:4000];

for iN=1:length(N)
    % Prolate spheroid
    bbechopdf(N(iN),mix_r,ns,gate_len,tx,sdir,'prosph_0.05',...
              [],2,0.05,[]);
end

