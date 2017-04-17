% 2012 11 08  Calculate echo pdf model


addpath /mnt/storage/prolatespheroid/
addpath /mnt/storage/modeling_code_current/

%{
% No beampattern
fname = 'tx_nosys';
tx_opt = 2;  % ideal tx, no system response
taper_opt = 0;
N = [10,300];
for iN=1:length(N)
    model_n_scat_no_bp(N(iN),sampleN,gate_len,tx_opt,...
                       taper_opt,save_opt,save_dir,fname);
end
%}

%{
% With beampattern
fname = 'tx_nosys';
tx_opt = 2;  % ideal tx, no system response
taper_opt = 0;
N = [1000];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gate_len,tx_opt,...
                            taper_opt,save_opt,save_dir,fname);
end
%}

%{
% Middle Hann window
fname = 'midhann_mid';
tx_opt = 1;  % square chirp
taper_opt = 4;  % mid Hann mid
N = [10,300];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gate_len,tx_opt,...
                            taper_opt,save_opt,save_dir,fname);
end

fname = 'midhann_narrow';
tx_opt = 1;  % square chirp
taper_opt = 7;  % mid Hann narrow
N = [10,300];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gate_len,tx_opt,...
                            taper_opt,save_opt,save_dir,fname);

end

fname = 'midhann_wide';
tx_opt = 1;  % square chirp
taper_opt = 9;  % mid Hann wide
N = [10,300];
for iN=1:length(N)
    model_n_scat_freqdep_bp(N(iN),sampleN,gate_len,tx_opt,...
                            taper_opt,save_opt,save_dir,fname);

end
%}

%{
fname = 'mixed';
tx_opt = 3;  % actual tx
taper_opt = 0;  % no additional taper
Nw = 1000;
Ns = [5,500];
r_sw = [1,10];
for iNS=1:length(Ns)
    N = [Nw,Ns(iNS)];
    model_n_scat_freqdep_bp_mixed(N,r_sw,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,0);
end

Nw = 1000;
Ns = [5,500];
r_sw = [1,5];
for iNS=1:length(Ns)
    N = [Nw,Ns(iNS)];
    model_n_scat_freqdep_bp_mixed(N,r_sw,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,0);
end

Nw = 1000;
Ns = [5,500];
r_sw = [1,30];
for iNS=1:length(Ns)
    N = [Nw,Ns(iNS)];
    model_n_scat_freqdep_bp_mixed(N,r_sw,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,0);
end
%}

%{
% monotype aggregation, with resp_x for noise-added models
gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121116_pdfmodel_w_respenv_5e4';
fname = 'tx_sys';
tx_opt = 3;  % actual tx
taper_opt = 0;  % no additional taper
%N = [10,50,300,500:100:1000,20:10:40,60:10:90,100,200,400];
N = [5:10:95,150:100:950];
r = 1;
for iN=1:length(N)
    model_n_scat_freqdep_bp_mixed(N(iN),r,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,0);
end
%}

% monotype aggregation, no system response
gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121127_tx_nosys_bbmodel';
fname = 'tx_nosys';
tx_opt = 2;  % ideal tx
taper_opt = 0;  % no additional taper
%N = [10,50,300,500:100:1000,20:10:40,60:10:90,100,200,400];
N = [10,300,50];
r = 1;
for iN=1:length(N)
    model_n_scat_freqdep_bp_mixed(N(iN),r,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,0);
end


% monotype aggregation, high-freq Hann
gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121127_hannwin_hf_lf';
fname = 'hann_hf';
tx_opt = 1;  % square chirp
taper_opt = 11;  % high-freq Hann, 20kHz bandwidth
%N = [10,50,300,500:100:1000,20:10:40,60:10:90,100,200,400];
N = [10,300,50];
r = 1;
for iN=1:length(N)
    model_n_scat_freqdep_bp_mixed(N(iN),r,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,0);
end

% monotype aggregation, low-freq Hann
gate_len = 0.5;
sampleN = 5e4;
save_opt = 1;
save_dir = '/mnt/storage/ECHO_STAT/20121127_hannwin_hf_lf';
fname = 'hann_lf';
tx_opt = 1;  % square chirp
taper_opt = 12;  % low-freq Hann, 20kHz bandwidth
%N = [10,50,300,500:100:1000,20:10:40,60:10:90,100,200,400];
N = [10,300,50];
r = 1;
for iN=1:length(N)
    model_n_scat_freqdep_bp_mixed(N(iN),r,sampleN,gate_len,...
                 tx_opt,taper_opt,save_opt,save_dir,fname,0);
end
