%SUMMARY
% Author: Connor Gallimore
% 12/28/2023

% This function creates a matrix of statistical values by giving the 
% 'imagesc()' function the appearance of a table using grid lines. 
% The 'directional' aspect of its name refers to the fact that it can 
% accept negative values, which, when they exceed the magnitude of 
% 'stat_thresh', can be plotted as a different color to offer more rich
% information about the statistical results.

% Required arguments: 
    % 'stat_matrix', a 2D m-by-n matrix spanning all real numbers
    % 'stat_thresh', the threshold at which more extreme values will be
    %                considered significant

% Optional Name,Value pairs
    % 'rowlabels', n-by-1 or 1-by-n string array.
    % 'collabels', n-by-1 or 1-by-n string array.
    % 'pvals',     a 2 element cell array where {1} contains 'p_matrix', a
    %              matrix of equivalent size as 'stat_matrix' with 
    %              corresponding p-values, and {2} is 'p_thresh' for
    %              significance.  This argument overlays the p-values as
    %              text, rather than the computed statistical value
    % 'XAxisLocation', followed by 'bottom', 'top', or 'origin'
    % 'YAxisLocation', followed by 'left', 'right', or 'origin'

% Outputs:
% a structure array of handles 'h' to axes, major/minor grid lines,
% row/col labels, table_image, text_labels overlaying the values, and 
% colorbar
%--------------------------------------------------------------------------


function h = statHeatmapDirectional(stat_matrix, stat_thresh, varargin)

tmp_cl=   strcmpi(varargin, 'collabels'); 
tmp_rl=   strcmpi(varargin, 'rowlabels'); 
tmp_pval= strcmpi(varargin, 'pvals'); 
tmp_xloc= strcmpi(varargin, 'XAxisLocation');
tmp_yloc= strcmpi(varargin, 'YAxisLocation');


% to use imagesc
xlen= size(stat_matrix, 2);
ylen= size(stat_matrix, 1);

if any(tmp_rl);   y_labels= varargin{find(tmp_rl) + 1};
else;             y_labels= strings(ylen, 1); % default row labels empty strings
end

if any(tmp_cl);   x_labels= varargin{find(tmp_cl) + 1};
else;             x_labels= strings(xlen, 1); % default col labels empty strings
end

stat_str= strings(ylen, xlen); 

if any(tmp_pval)
    p_vars= varargin{find(tmp_pval) + 1}; 
    p_matrix= p_vars{1};
    p_thresh= p_vars{2};
    p_criteria= p_matrix < p_thresh;
    stat_str(stat_matrix < -stat_thresh & p_criteria)= string(p_matrix(stat_matrix < -stat_thresh & p_criteria)); 
    stat_str(stat_matrix > stat_thresh & p_criteria)= string(p_matrix(stat_matrix > stat_thresh & p_criteria)); 
else
    stat_str(stat_matrix < -stat_thresh)= string(abs(stat_matrix(stat_matrix < -stat_thresh))); 
    stat_str(stat_matrix > stat_thresh)=  string(stat_matrix(stat_matrix > stat_thresh)); 
end

[yt, xt]= ndgrid(1:ylen, 1:xlen); 

% plot data and x/y tick labels
stat_image= imagesc(stat_matrix);
h.axes= gca;
set(h.axes, 'XTick', 1:length(x_labels), 'YTick', 1:length(y_labels)); 
h.axes.TickLength(1)= 0;
xticklabels(x_labels);
yticklabels(y_labels); 

cb= colorbar(h.axes);  

hold on 

% plot major grid lines
h.major_xlines= arrayfun(@(x) xline(h.axes, x, 'k-', 'Alpha', 1), 0.5:1:(xlen + 0.5));
h.major_ylines= arrayfun(@(y) yline(h.axes, y, 'k-', 'Alpha', 1), 0.5:1:(ylen + 0.5));

axis image

% overlay text
th= text(xt(:), yt(:), stat_str(:), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'Center');

if any(tmp_xloc)
    set(h.axes, varargin{tmp_xloc}, varargin{find(tmp_xloc) + 1}); 
end

if any(tmp_yloc)
    set(h.axes, varargin{tmp_yloc}, varargin{find(tmp_yloc) + 1}); 
end

h.table_image= stat_image; 
h.text_labels= th;
h.colorbar= cb;

end