function [x,y,varargout] = findEchoDist(data,varargin)
%function [xlog,nlog,xlin,nlin,hit_num] = findEchoDist(data,opt)

% Get normalized amplitude distribution
% INPUT
%   data        echo amp data
%   varargin    varargin{1}: #points on the x-axis; default=100
%               varargin{2}: 0-log, 1-lin; default=0
%               varargin{3}: 1-plot, 0-no plot; default=0
%
% OUTPUT
%   x           x-axis
%   y           pdf
%   varargout   # of hits in each x-axis bin
%
% WJL 2010/05/16
% WJL 2010/07/19 modified
% WJL 2010/10/18 modified take out the normalization for raw distribution
% WJL 2012/01/02 modified to use varagrout
% WJL 2012/04/05 modified to use varargin
% WJL 2012/11/01 simplify the code...

if ~isempty(varargin)
    if length(varargin)==3
        npt = varargin{1};
        ll_opt = varargin{2};
        plot_opt = varargin{3};
    elseif length(varargin)==2
        npt = varargin{1};
        ll_opt = varargin{2};
        plot_opt = 0;
    else
        npt = varargin{1};
        ll_opt = 0;
        plot_opt = 0;
    end
else
    npt = 100;
    plot_opt = 0;
end

if ll_opt==0
    edges = logspace(log10(min(data)),log10(max(data)),npt+1);
    x = sqrt(edges(1:end-1).*edges(2:end));
else
    edges = linspace(min(data),max(data),npt+1);
    x = mean([edges(1:end-1);edges(2:end)],1);
end

hitnum = histc(data,edges);
y = hitnum(1:end-1)./diff(edges);
y = y/trapz(x,y);

x = x.';
y = y.';
varargout{1} = hitnum;

if plot_opt==1
    rayl = raylpdf(x,1/sqrt(2));
    figure;
    hr = loglog(x,rayl,'k');
    hold on
    hy = loglog(x,y,'rx');
    legend({'Rayleigh','Data'});
    xlabel('Echo amplitude','fontsize',14);
    ylabel('PDF','fontsize',14);

    ylim([1e-5 1e2])
    xlim([1e-3 1e2])
end




