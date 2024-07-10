% SUMMARY
% Author: Nils Haentjens
% Created: Jan 15, 2018
% Modified: Connor Gallimore, Dec 5, 2022

% Scatter plot where each point is colored by the spatial density of nearby
% points. The function use the kernel smoothing function to compute the
% probability density estimate (PDE) for each point. It uses the PDE has
% color for each point.
%
% Input
%     x <Nx1 double> position of markers on X axis
%     y <Nx1 double> posiiton of markers on Y axis
%     varargin can be used to send a set of instructions to the scatter function
%           Supports the MarkerSize parameter
%           Does not support the MarkerColor parameter

%     'DenseSort' can be specified to ensure the highest probability values
%     are plotted on top
%
% Output:
%     h returns handles to the scatter objects created
%
% Example
%     % Generate data
%     x = normrnd(10,1,1000,1);
%     y = x*3 + normrnd(10,1,1000,1);
%     % Plot data using probability density estimate as function
%     figure(1); 
%     scatter_kde(x, y, 'filled', 'MarkerSize', 100);

%     alternatively:
%     scatter_kde(x, y, 'filled', 'DenseSort');

%     % Add Color bar
%     cb = colorbar();
%     cb.Label.String = 'Probability density estimate';
%--------------------------------------------------------------------------


function varargout = scatter_kde(x, y, varargin)

if isrow(x)
    x = x(:);
end
if isrow(y)
    y = y(:);
end

% Use Kernel smoothing function to get the probability density estimate (c)

% REMEMBER that a PDF is *NOT* a probability, and values here can exceed 1
% e.g. "the uniform distribution on the interval [0, ½] 
% has probability density f(x) = 2 for 0 ≤ x ≤ ½ and f(x) = 0 elsewhere."
c = ksdensity([x,y], [x,y]);

if nargin > 2
  % Set Marker Size
  mSz = strcmpi(varargin, 'MarkerSize');
  if any(mSz)
      MarkerSize = varargin{find(mSz)+1}; 
      varargin(mSz:mSz+1) = [];
  else
      MarkerSize = [];
  end
  % Sort so that the highest probabilities are plotted on top
  dSort = strcmpi(varargin, 'DenseSort');
  if any(dSort)
      [c, ix]= sort(c); 
      x= x(ix);
      y= y(ix); 
      varargin(dSort) = [];
  end

  vals_only= strcmpi(varargin, 'ValsOnly'); 
  if any(vals_only)
      h= []; 
  else
      % Plot scatter plot
      h = scatter(x, y, MarkerSize, c, varargin{:});
  end
  
else
  h = scatter(x, y, [], c);
end

varargout= { h, x, y, c }; 

end