% 2011 12 08
% Check again on the ir of beampattern response
% because the previous one had a constant error in it
% Need to make sure everything is correct

df = 1e-4;
const = 3;
f = 0:df:100;
bp = 2*besselj(1,const*2*pi*f)./(const*2*pi*f);
bp(1) = 1;

dt = 1/(2*max(f));
t = -const:dt:const;
bp_t = 2/abs(const)*sqrt(2/pi)*sqrt(1-(t/const).^2);

bp_full = [bp,fliplr(conj(bp(2:end)))];
bp_ifft = length(bp_full)*ifftshift(ifft(bp_full))*df;
bp_ifft = bp_ifft*sqrt(2*pi);  % sqrt(2*pi) is for the conversion between
                               % differet forms of fourier transform
t_ifft = ((1:length(bp_ifft))-length(bp))*dt;
bpt_half_len = floor((length(bp_full)+1)/2);

figure;
plot(t_ifft,abs(bp_ifft));
hold on
plot(t,bp_t,'r--');
xlim([-5 5])

%% use ifft to obtain ir
%bf_both = [bf;conj(flipud(bf(2:end,:)))];
%bt_ifft = length(bf_both)*ifft(bf_both);
%bt_ifft = fftshift(bt_ifft)*df;
%bL = floor(size(bf_both,1)/2);

%dt_bp = 1/(2*max(freq));
%t_bp = dt_bp*((0:size(bt_ifft,1)-1)-bL);
%[~,ind1] = min(abs(2*const+t_bp));
%[~,ind2] = min(abs(2*const-t_bp));
%t = t_bp(ind1-1:ind2+1); 
%bt = bt_ifft(ind1-1:ind2+1); 
