% SUMMARY
% Author: Connor Gallimore
% 12/06/2022

% Constructs a bivariate probability heatmap for the relationship between
% two variables. Values can be computed in terms of counts or probability

% Parameters as 'name',value pairs:
% - 'nbins': Array in the form of [nxbins nybins] to set the
% number of bins in each dimension
% - 'edges': Cell in the form of {[x_edges] [y_edges]} to set
% custom bin edges for each dimension
%
% For more details, refer to hist3

% Custom:
% - 'probability': converts the colormap from counts to probability

% Modified from: densityplot
% https://www.mathworks.com/matlabcentral/fileexchange/65166-densityplot-x-y-varargin
%--------------------------------------------------------------------------


function varargout = densityContour(x, y, varargin)

% Process x and y inputs
if isrow(x)
    x = x(:);
end
if isrow(y)
    y = y(:);
end
X= [x, y];   % combine for hist3 input

% Process params
arg_bins= strcmpi(varargin, 'nbins'); 
arg_edge= strcmpi(varargin, 'edges'); 
arg_ncnt= cellfun(@isscalar, varargin);  % n contours, controls smoothness
arg_norm= strcmpi(varargin, 'probability');
arg_vals= strcmpi(varargin, 'valsOnly'); 

if any(arg_bins)
    locB= find(arg_bins); 
    hist_args= {varargin{locB}, varargin{locB+1}}; 
else
    hist_args= {};
end

if any(arg_edge)
    locE= find(arg_edge);
    hist_args= horzcat(hist_args, varargin{locE}, varargin{locE+1}); 
end

% get bivariate histogram
[N, bin_ctrs]= hist3(X, hist_args{:}); 

if any(arg_norm)     % returns probability rather than counts, if specified
    N= N ./ sum(N(:)); % such that sum of all bar heights/magnitudes = 1
end

% set number of contours
if any(arg_ncnt)
    n_contours= varargin{arg_ncnt}; 
else
    n_contours= min(size(N));   % defaults to lowest bin dim size
end

if any(arg_vals)
    h= []; M= []; 
else
    [M, h]= contourf(bin_ctrs{1}, bin_ctrs{2}, N', n_contours, 'linecolor', 'none'); 
end

varargout= {h, M, N, bin_ctrs}; 

end

