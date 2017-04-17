function win = win_chirp(taper_opt,y)

% windowing the chirp
if taper_opt==1  % full Gaussian taper
    winL = round(length(y)/2);
    win = gausswin(winL);
    win = [win(1:floor(winL/2))',ones(1,length(y)-winL),win(floor(winL/ ...
                                                      2)+1:end)'];
elseif taper_opt==2  % HF Hann taper
    winL = round(length(y)/2);
    win = hann(winL);
    win = [zeros(1,length(y)/2),win(1:floor(winL/2))',ones(1,winL/2)];
elseif taper_opt==3  % LF Hann taper
    winL = round(length(y)/2);
    win = hann(winL);
    win = [ones(1,winL/2),win(floor(winL/2)+1:end)',zeros(1,length(y)/2)];
elseif taper_opt==4  % middle Hann taper
    winL = round(length(y)/4);
    win = hann(winL);
    win = [zeros(1,winL),win(1:winL/2)',ones(1,winL),win(winL/2+1:end)',zeros(1,winL)];
elseif taper_opt==5  % upper Hann taper
    winL = round(length(y)/4);
    win = hann(winL);
    win = [zeros(1,winL*2),win(1:winL/2)',ones(1,winL),win(winL/2+1:end)'];
elseif taper_opt==6  % lower Hann taper
    winL = round(length(y)/4);
    win = hann(winL);
    win = [win(1:winL/2)',ones(1,winL),win(winL/2+1:end)',zeros(1,winL*2)];
elseif taper_opt==7  % middle Hann taper narrow1
    winL = round(length(y)/4);
    win = hann(winL);
    win = [zeros(1,winL*11/8),win(1:winL/2)',ones(1,winL/4),win(winL/2+1:end)',zeros(1,winL*11/8)];
elseif taper_opt==77  % middle Hann taper narrow1 very narrow  1/4 y len
    winL = round(length(y)/4);
    win = hann(winL);
    win = [zeros(1,winL*3/2),win',zeros(1,winL*3/2)];
elseif taper_opt==777  % middle Hann taper narrow1 very narrow  1/16 y len
    winL = round(length(y)/16);
    win = hann(winL);
    win = [zeros(1,winL*15/2),win',zeros(1,winL*15/2)];
elseif taper_opt==8  % middle Hann taper narrow2 --> this one not used
    winL = round(length(y)/4);
    win = hann(winL);
    win = [zeros(1,winL*5/4),win(1:winL/2)',ones(1,winL/2),win(winL/2+1:end)',zeros(1,winL*5/4)];
elseif taper_opt==9  % middle Hann taper wider
    winL = round(length(y)/4);
    win = hann(winL);
    win = [zeros(1,winL/2),win(1:winL/2)',ones(1,winL*2),win(winL/2+1:end)',zeros(1,winL/2)];
elseif taper_opt==10  % middle Hann taper wider2
    winL = round(length(y)/4);
    win = hann(winL);
    win = [win(1:winL/2)',ones(1,winL*3),win(winL/2+1:end)'];
elseif taper_opt==11  % upper Hann taper 20kHz band
    winL = round(length(y)/4);
    win = hann(winL);
    win = [zeros(1,winL),win(1:winL/2)',ones(1,winL*2),win(winL/2+1:end)'];
elseif taper_opt==12  % lower Hann taper 20kHz band
    winL = round(length(y)/4);
    win = hann(winL);
    win = [win(1:winL/2)',ones(1,winL*2),win(winL/2+1:end)',zeros(1, ...
                                                      winL)];
else
    win = ones(size(y));
end
