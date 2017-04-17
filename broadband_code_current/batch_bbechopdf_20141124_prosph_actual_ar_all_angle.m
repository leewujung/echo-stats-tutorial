% 2014 04 27  Figure out where the code was left in August 2013
% 2014 05 26  Run more fish scatterer models for comparison
% 2014 09 04  Run bbechopdf using prolate spheroid scatterers with
%             realistic aspect ratio
% 2014 11 23  Found out an error in prosph_amp.m which affects all
%             results calculated before using all angles of
%             orientation. But the models with ar=0.05 currently
%             used in the paper were calculated using prescribed
%             angle of orientation distribution so are not affected
%             bu this error.
% 2014 11 25  Use [0,2*pi] uniform as the default angle of orientation


mix_r = 1;
ns = 5e4;
gate_len = 0.5;
save_opt = 1;
tx = 3;
sdir = '/mnt/storage/ECHO_STAT/20141125_bbechopdf_prolate_sphroid_realistic_ar_single_length_all_angle';
if ~exist(sdir,'dir')
    mkdir(sdir);
end

N=[20,100,1000,1500,2000,2500,3000]';
%N = [10:10:90,100:100:4000];
%N = [10,50,300];

for iN=1:length(N)
    % Prolate spheroid
    bbechopdf_20141124(N(iN),mix_r,ns,gate_len,tx,sdir,'prosph_0.05',...
              [],2,[0.05,0.25,1],[NaN,NaN]);
end


N = [40,60:10:90,200,400:100:900,1100:100:1400,1600:100:1900,2100:100:2400,2600:100:4000];
%N = [10:10:90,100:100:4000];
%N = [10,50,300];

for iN=1:length(N)
    % Prolate spheroid
    bbechopdf_20141124(N(iN),mix_r,ns,gate_len,tx,sdir,'prosph_0.05',...
              [],2,[0.05,0.25,1],[NaN,NaN]);
end

