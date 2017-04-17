function [y,t_y] = gen_tx(tx_opt,varargin)
% tx_opt   1-square linear chirp
%          2-ideal transmit signal (no system response)
%          3-actual transmit signal (with system response)

if tx_opt==1 & nargin~=1
    disp('Invalid inputs!');
    return;
elseif tx_opt==1 & nargin==1   % Square chirp signal
    t_y = 0:3.6e-6:0.01;    % pulse length [sec]
    y = chirp(t_y,30e3,0.01,70e3);   % start/end freq [Hz]
elseif tx_opt==2  % Ideal transmit signal (no system resp)
    [y,t_y] = chirp_no_sys(120,100);
    y = y/max(y);
elseif tx_opt==3  % Actual transmit signal (with system resp)
    [y,t_y] = chirp_w_sys(120,100);
    y = y/max(y);
end
