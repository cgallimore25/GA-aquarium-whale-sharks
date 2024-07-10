%SUMMARY
% Author: Connor Gallimore
% 01/03/2024

% This function creates a dot plot representing a matrix of statistical
% values, comparable to R's 'corrplot' function for many pairwise/cross
% correlations. More generally, it can be used for any matrix containing
% results like p-values or statistical tests.

% inspired by:
% https://www.mathworks.com/matlabcentral/answers/699755-fancy-correlation-plots-in-matlab
% and
% http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2

% Required arguments: 
    % 'matrix', a 2D m-by-n matrix

% Optional Name,Value pairs
    % 'lower', 'upper', or 'full', specifying the plot type
    %              for lower/upper, acceptable values to follow are 0, 1, 
    %              or -1 to keep/remove upper/lower identity diagonals, 
    %              respectively.
    % 'majorgrid', 1 or 0, toggles black lines encasing all points in the
    %              matrix.
    % 'minorgrid', 1 or 0, toggles light gray lines intersecting data pts.
    % 'dotscale',  a scalar controlling the size of plotted dots.
    % 'rowlabels', n-by-1 or 1-by-n string array.
    % 'collabels', n-by-1 or 1-by-n string array.
    % 'alpha',     numeric scalar between the closed interval [0-1]
    %              indicating opacity of the plotted dots. 
    % 'overlayvals', followed by 'all' or 'sigonly', telling the function
    %              whether to display the values in the matrix. 'sigonly'
    %              must be followed up with a second numeric scalar
    %              representing the threshold, whereby values below that
    %              threshold will not have their text overlaid. 
    % 'precision', integer specifying rounding precision for text labels 
    %              (i.e. the max number of decimal places).

% Outputs:
% a structure array of handles 'h' to axes, dots, major/minor grid lines,
% row/col labels, text overlaying the values, and colorbar

% Examples: 
% see 'xcorrDotPlot_demo.mlx' for usage tips and tricks
%--------------------------------------------------------------------------


function h = xcorrDotPlot(matrix, varargin)

sz= size(matrix);
[~, loc]= max(sz);

% make general coordinate system
if sz(1) == sz(2)
    n= sz(1);
    y_coords= repmat(n+1, n, n) - (1:n)';
    x_coords= repmat(1:n, n, 1);
else
    [y_coords, x_coords]= ndgrid(sz(1):-1:1, 1:sz(2));
end

% variable args
tmp_tri= strcmpi(varargin, 'lower') | strcmpi(varargin, 'upper') | strcmpi(varargin, 'full');
tmp_gr=  strcmpi(varargin, 'majorgrid'); 
tmp_mg=  strcmpi(varargin, 'minorgrid'); 
tmp_fa=  strcmpi(varargin, 'alpha'); 
tmp_rl=  strcmpi(varargin, 'rowlabels');
tmp_cl=  strcmpi(varargin, 'collabels');
tmp_do=  strcmpi(varargin, 'dotscale'); 
tmp_va=  strcmpi(varargin, 'overlayvals'); 
tmp_pr=  strcmpi(varargin, 'precision'); 

if any(tmp_tri);  plottype= lower(varargin{tmp_tri}); 
else;             plottype= 'full';
end

if ~strcmpi(plottype, 'full');  k= varargin{find(tmp_tri) + 1};
else;                           k= 0;   % default keeps identity diagonal
end

if any(tmp_gr);   g= varargin{find(tmp_gr) + 1};
else;             g= 1;     % major grid ON by default
end

if any(tmp_mg);   mg= varargin{find(tmp_mg) + 1};
else;             mg= 0;    % minor grid OFF by default
end

if any(tmp_fa);   a= varargin{find(tmp_fa) + 1};
else;             a= 0.8;   % default marker alpha
end

if any(tmp_rl);   RL= varargin{find(tmp_rl) + 1};
else;             RL= strings(sz(1), 1); % default row labels empty strings
end

if any(tmp_cl);   CL= varargin{find(tmp_cl) + 1};
else;             CL= strings(sz(2), 1); % default col labels empty strings
end

if any(tmp_do);   dotSc= varargin{find(tmp_do) + 1};
else;             dotSc= 1; % default dot scaling is 1 (i.e. none)
end

if any(tmp_pr);   prec= varargin{find(tmp_pr) + 1};
else;             prec= 2; 
end

% pre-compute all coordinates for lines and data points
if sz(1) == sz(2)
    switch plottype
        case 'lower'
    
            data= tril(matrix, k);
            data(~tril(matrix, k))= nan; 
            
            % x and y coords
            yc= tril(y_coords, k) + 0.5;
            xc= tril(x_coords, k) + 0.5;
            xc(xc == 0.5) = NaN;
            
            [xtr, ytr, xtc, ytc, maj_gr, min_gr, RL, CL]= adjustLineAndLabelCoords(plottype, k, n, RL, CL);
    
        case 'upper' 
    
            data= triu(matrix, k);
            data(~triu(matrix, k))= nan; 
    
            yc= triu(y_coords, k) + 0.5;
            xc= triu(x_coords, k) + 0.5;
            xc(xc == 0.5) = NaN;
    
            [xtr, ytr, xtc, ytc, maj_gr, min_gr, RL, CL]= adjustLineAndLabelCoords(plottype, k, n, RL, CL);

        otherwise
            data= matrix;
            data(~matrix)= nan;
            yc= y_coords + 0.5; 
            xc= x_coords + 0.5;
            xtr= ones(1, n) * 0.5; 
            ytr= (n:-1:1) + 0.5;
            xtc= (1:n) + 0.5;
            ytc= repmat(n+1, 1, n);
    end
else
    data= matrix;
    data(~matrix)= nan;
    yc= y_coords + 0.5; 
    xc= x_coords + 0.5;
    xtr= ones(1, sz(1)) * 0.5; 
    ytr= (sz(1):-1:1) + 0.5;
    xtc= (1:sz(2)) + 0.5;
    ytc= repmat(sz(2)+1, 1, sz(2));
end

% axis handle
h.axes= gca; hold on

% plot minor grid lines
if logical(mg) && ~strcmpi(plottype, 'full')
    h.minor_ylines= line(min_gr(1:2, :), min_gr(3:4, :), 'color', [.8 .8 .8]);   % horizontal minor lines
    h.minor_xlines= line(min_gr(3:4, :), min_gr(1:2, :), 'color', [.8 .8 .8]);   % vertical minor lines
elseif logical(mg) && strcmpi(plottype, 'full')
    mxh= [xc(:, 1)-0.5, xc(:, end)+0.5]';
    myh= [yc(:, 1) yc(:, end)]';
    h.minor_ylines= line(mxh, myh, 'color', [.8 .8 .8]);   % horizontal minor lines
    mxv= xc(1:2, :);
    myv= [ones(1, sz(2)); repmat(sz(1)+1, 1, sz(2))];
    h.minor_xlines= line(mxv, myv, 'color', [.8 .8 .8]);   % vertical minor lines
end

% plot data points
h.points= scatter(xc(:), yc(:), abs(data(:)) * dotSc, data(:), 'filled', 'MarkerFaceAlpha', a);

% plot major grid lines
if logical(g) && ~strcmpi(plottype, 'full')
    h.major_ylines= line(maj_gr(1:2, :), maj_gr(3:4, :), 'color', 'k');   % horizontal major lines
    h.major_xlines= line(maj_gr(3:4, :), maj_gr(1:2, :), 'color', 'k');   % vertical major lines
elseif logical(g) && strcmpi(plottype, 'full')
    h.major_ylines= arrayfun(@(y) yline(h.axes, y, 'k-', 'Alpha', 1), 1:sz(1)+1);
    h.major_xlines= arrayfun(@(x) xline(h.axes, x, 'k-', 'Alpha', 1), 1:sz(2)+1);
end

% add row/col labels
h.row_lbls= text(xtr, ytr, RL, 'HorizontalAlignment', 'right');
h.col_lbls= text(xtc, ytc, CL, 'HorizontalAlignment', 'right', 'Rotation', 270);

% overlay value text
if any(tmp_va)
    val_type= lower(varargin{find(tmp_va) + 1});
    val_txt=  string(round(data, prec));
    switch val_type
        case 'all'
            h.val_txt= text(xc(:), yc(:), val_txt(:), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'Center');
        case 'sigonly'
            sig_thresh= abs(varargin{find(tmp_va) + 2}); % make thresh (+)
            criteria=   data <= sig_thresh & data >= -sig_thresh;
            val_txt(criteria)= ""; 
            h.val_txt= text(xc(:), yc(:), val_txt(:), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'Center');
    end
end

h.cbar= colorbar(h.axes);
h.axes.Visible = 'off';
h.axes.Position(4) = h.axes.Position(4)*((sz(loc)-1)/sz(loc));
axis(h.axes, 'image')
set(gcf, 'color', 'w');


end


%% HELPER FUNCTIONS--------------------------------------------------------

function [xtr, ytr, xtc, ytc, major_g, minor_g, rowLabels, colLabels]= ...
          adjustLineAndLabelCoords(triangle_type, k, n, rowLabels, colLabels)

switch triangle_type
    case 'lower'
        xl = [2:n+1; ones(1, n)];
        
        if k <= 0
            % major x and y grid lines
            xl = [xl(:, 1:end+k), xl(:, end+k)];
            yl = repmat((n+1)+k:-1:1, 2, 1);

            % x and y row/col label coords
            xtr= ones(1, n+k) * 0.5;
            ytr= (n+k:-1:1) + 0.5;
            xtc= (2:(n+1)+k) - 0.5; 
            ytc= ((n+1)+k:-1:2) + 0.5; 

            % row/col labels
            rowLabels= rowLabels(1-k:end); 
            colLabels= colLabels(1:end+k); 

        elseif k > 0
            xl = [xl, xl(:, end)]; 
            xl(1, 1:n-k) = xl(1, 1:n-k) + k;
            xl(1, n-(1:k-1))= xl(1, n-(1:k-1)) + (1:(k-1));
            yl = repmat((n+1):-1:1, 2, 1); 
            
            % x and y row/col label coords
            xtr= ones(1, n) * 0.5;
            ytr= (n:-1:1) + 0.5;
            xtc= (2:(n+1)) - 0.5; 
            ytc= [repmat((n+1), 1, k+1), n:-1:2+k] + 0.5;

        end

        % minor x and y grid lines
        mx= xl(:, 1:end-1);
        my= yl(:, 1:end-1) - 0.5;

    case 'upper'
        if k >= 0
            xl = [1+k:n+1; repmat(n+1, 1, (n+1)-k)];
            xl = [xl(:, 1), xl(:, 1:end-1)];
            yl = repmat(n+1:-1:1+k, 2, 1);
    
            xtr= 1+k:n; 
            ytr= (n:-1:1+k) + 0.5; 
            xtc= (1+k:n) + 0.5;
            ytc= repmat(n+1, 1, n-k);

            rowLabels= rowLabels(1:end-k); 
            colLabels= colLabels(1+k:end); 
        elseif k < 0
            xl = [1:n+1; repmat(n+1, 1, (n+1))];
            xl = [repmat(xl(:, 1), 1, -k+1), xl(:, 1:end-1+k)];
            yl = repmat(n+1:-1:1, 2, 1); 

            xtr= [ones(1, 1-k), 2:n+k];
            ytr= (n:-1:1) + 0.5; 
            xtc= (1:n) + 0.5;
            ytc= repmat(n+1, 1, n);
        end

        mx= xl(:, 2:end);
        my= yl(:, 2:end) + 0.5;
end

major_g=    [xl; yl];
minor_g=    [mx; my];

end
