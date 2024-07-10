% SUMMARY
% Author: Connor Gallimore
% 10/09/2023

%   Add summary statistic lines to output of 'plotSpread... .m' functions 
%   to display mean or median over the scatter distributions

% Optional arguments (varargin)
%   specify line color using an rgb array where each element goes from 0-1
%   default is black ('k'; or equivalently [0 0 0]) 

%   specify width of the line using a scalar, which will be distributed
%   evenly on either side.  default is 0.5

%   specify line thickness using a scalar.  default is 2

% Optional outputs (varargout)
%   1 output argument will supply the handle to the overlaid lines
%   2 output arguments will supply values of the summary statistic
%--------------------------------------------------------------------------


function varargout = addSummaryStatLines(data, ngroups, summaryStat, varargin)

% process optional Name,Val pairs, assign defaults if absent
tmpC= strcmpi(varargin, 'color');
tmpW= strcmpi(varargin, 'width'); 
tmpL= strcmpi(varargin, 'linewidth'); 
tmpX= strcmpi(varargin, 'xVals');

if any(tmpC);  color= varargin{find(tmpC) + 1}; 
else;          color= 'k'; 
end

if any(tmpW);  width= varargin{find(tmpW) + 1}; 
else;          width= 0.5; 
end

if any(tmpL);  linew= varargin{find(tmpL) + 1}; 
else;          linew= 2; 
end

if any(tmpX);  xt= varargin{find(tmpX) + 1}; 
               xc= [xt + width/2; xt - width/2]; 
else;          xc= [(1:ngroups) + width/2; (1:ngroups) - width/2]; 
end

% parse if matrix or cell input? 

if iscell(data)
    % calculate chosen summary statistic for cell
    if strcmpi(summaryStat, 'mean')
        vals= cellfun(@(x) mean(x(:), 'omitnan'), data); 
    elseif strcmpi(summaryStat, 'median')
        vals= cellfun(@(x) median(x(:), 'omitnan'), data);
    end
else
    % get dim to calculate along
    dim= find( ~(size(data) == ngroups) );
    
    % calculate chosen summary statistic for matrix
    if strcmpi(summaryStat, 'mean')
        vals= mean(data, dim, "omitnan");
    elseif strcmpi(summaryStat, 'median')
        vals= median(data, dim, "omitnan");
    end
end

yc= [vals; vals]; 

statLines= plot(xc, yc, '-', 'Color', color, 'LineWidth', linew);

varargout = {statLines, vals};

end