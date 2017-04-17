

ORI_D = '/mnt/storage/broadband_echopdf_models/fish_defaultLenDistr_mean5_std20/';
NEW_D = '/mnt/storage/broadband_echopdf_models/fish_defaultLenDistr_mean5_std20_new/';
files = dir([ORI_D,'/fish_*.mat']);

for iF=1:length(files)
    fname = files(iF).name;
    fname_new = [fname(1:30),fname(32:end)];
    copyfile([ORI_D,'/',fname],[NEW_D,'/',fname_new]);
end

