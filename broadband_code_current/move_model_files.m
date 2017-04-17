% 2012 11 08  Move and rename model files

ORI_D = '/mnt/storage/ECHO_STAT/20120731_bp_results';
NEW_D = '/mnt/storage/ECHO_STAT/20121108_bbfig_models';
files = dir([ORI_D,'/middle_hann_taper_narrow*.mat']);

newh = 'midhann_narrow';
for iF=1:length(files)
    [~,rem] = strtok(files(iF).name,'_');
    [~,rem] = strtok(rem,'_');
    [~,rem] = strtok(rem,'_');
    [~,rem] = strtok(rem,'_');
    nf = [newh,rem];
    copyfile([ORI_D,'/',files(iF).name],[NEW_D,'/',nf]);
end

