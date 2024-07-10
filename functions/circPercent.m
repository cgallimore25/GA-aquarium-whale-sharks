%SUMMARY
% Author: Connor Gallimore
% 12/26/2023

% This function creates arcs of a circle that represent percentages as the
% proportion of total circumference. I first saw this data visualization
% idea in a Nature paper (below) and wanted to implement it here in MATLAB.

% (Fig 4, Soula et al. 2023, https://www.nature.com/articles/s41593-023-01270-2)

% Required arguments:
    % 'data', a vector or matrix of proportions/percentages 
    % 'dim',  the dimension containing each component of the total

% Optional Name,Value pairs:
    % 'color', an m x 3 vector or matrix specifying RGB triplet(s)
    % 'orientation', 'horizontal' or 'vertical', for more than 1 data 
    %                series, determines whether they are plotted from left-
    %                to-right or top-to-bottom
    % 'precision', specifies the rounding precision for text labels (i.e. 
    %              the max number of decimal places)
    % 'radius', specifies the radius of the circle, indirectly impacting
    %           the width of the line arcs
    % 'scheme', 'category' or 'series', determines the coloring scheme. In
    %           one case, your 'color' matrix may represent the color you 
    %           want each common 'category' to be for all series (default).
    %           In another case, you may be specifying the base color you
    %           want your percentage components to be for each 'series'. In
    %           this option, subsequent percentages will be plotted darker


% Outputs:
% a structure of handles 'H' to the line 'arcs', text labels ('lbls'), and
% 'colors'

% Examples: 
% see 'circPercent_demo.mlx' for usage tips and tricks

% Additional notes:
% The origin and radius are arbitrary, as the width of line arcs will be 
% scaled as a function of radius and number of plot series, with axes 
% removed at the end anyways. I chose (0, 0) and 3, respectively.

% I used the helper function 'getAlignmentFromAngle' from MATLAB's 'pie.m'
% to position text labels, for which I take no credit and all rights go to
% Clay M. Thompson 3-3-94 and The MathWorks, Inc., Copyright 1984-2022 

% This function also uses 'distinguishable_colors.m' by Timothy E. Holy, 
% for which all rights are reserved to them (Copyright (c) 2010-2011, 
% https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors)
%--------------------------------------------------------------------------


function H = circPercent(data, dim, varargin)

% defaults -- center origin (h, k) and radius (r)
h= 0;
k= 0;
r= 3;

% parse Name,Val pairs
tmp_col= strcmpi(varargin, 'color'); 
tmp_ori= strcmpi(varargin, 'orientation');
tmp_pre= strcmpi(varargin, 'precision'); 
tmp_rad= strcmpi(varargin, 'radius'); 
tmp_sch= strcmpi(varargin, 'scheme'); 


if any(tmp_ori);  ori= varargin{find(tmp_ori) + 1};
else;             ori= 'horizontal';   % default horizontal
end

if any(tmp_sch);  sc= varargin{find(tmp_sch) + 1};
else;             sc= 'category';      % default color by category
end

if any(tmp_pre);  prec= varargin{find(tmp_pre) + 1};
else;             prec= 2;             % default round to second decimal
end

if any(tmp_rad);  r= varargin{find(tmp_rad) + 1};
end                          % overwrite default radius if specified

% anonymous functions validating proportions & optional color matrix input
validProps= @(x) isnumeric(x) && ismatrix(x) && all(all(x >= 0) & all(x <= 1));
validColor= @(x) validProps(x) && size(x, 2) == 3; 

% if groups distributed along columns, transpose
if dim ~= 2 
    data= data';
    dim= 2; 
end

n_dims= length(size(data)); 
nc_dim= dim; 
np_dim= find(~ismember(1:n_dims, dim));

np= size(data, np_dim);   % num percents to plot
nc= size(data, nc_dim);   % num components to total

% normalize groups that don't sum to 1
if any(sum(data, dim) > 1+sqrt(eps))
    ix2n= sum(data, dim) > 1+sqrt(eps);  % indices to normalize
    data(ix2n, :)= data(ix2n, :) ./ sum(data(ix2n, :), dim);
end

% make txt strings for later based on rounding precision
txt= strcat(string(round(data .* 100, prec)), repmat("%", np, nc)); 

% parse color / color scheme args
if any(tmp_col)
    col= varargin{find(tmp_col) + 1};
    if validColor(col)
        fc= col; 
    elseif isempty(col)
        if strcmpi(sc, 'series')
            fc= distinguishable_colors(np*2);
            fc= fc(np+1:end, :);
        else
            fc= distinguishable_colors(nc*2);
            fc= fc(nc+1:end, :);
        end
    else
        error('make sure your color matrix is either mx3 or left empty []')
    end
else
    fc= distinguishable_colors(nc*2);
    fc= fc(nc+1:end, :);
end
dark= linspace(1, 0, nc+1);  % darkness scaling
dark= dark(1:nc)'; 

% anonymous function for color decision rule
colorCat= @(x, y, z) strcmpi(x, 'category') && size(y, 1) == z; 

% resolve potential color / scheme input discrepancies
if ~colorCat(sc, fc, nc)
    if strcmpi(sc, 'series') && size(fc, 1) ~= np
        warning(['series scheme was indicated, but a color matrix was either unspecified or ' ...
                 'dim 1 of color matrix didnt match the number of data series(' num2str(np) ').' ...
                 ' making new colors as default.'])
        fc= distinguishable_colors(np*2);
        fc= fc(np+1:end, :);
        
    elseif strcmpi(sc, 'category') 
        warning(['category scheme was indicated, but dim 1 of color matrix didnt match the number ' ...
                 'of categories(' num2str(nc) '). making new colors as default.'])
        fc= distinguishable_colors(nc*2);
        fc= fc(nc+1:end, :);
    end
end


if np > 1  % if plotting a series of percentages, scale radius differently
    scF= 1.25;
else
    scF= 1.4; 
end

% adjust (h, k) for the desired orientation
if strcmpi(ori, 'horizontal')
    h= h:(3*r+1):(3*r+1)*np-1;      % step right for horizontal series
    k= zeros(1, np); 
else
    h= zeros(1, np); 
    k= k:-(3*r+1):-((3*r+1)*np-1);  % step down for vertical series
end

% pre-set ax limits
x1= h-r-np;   x2= h+r+np;    y1= k-r-np;    y2= k+r+np;
x_lo= x1(1);  x_hi= x2(np);  y_lo= y1(np);  y_hi= y2(1); 


% make arc start and finish points
a0= zeros(size(data, np_dim), 1);   % arbitrary starting point (0 degrees)
p360= data .* 360;
af= deg2rad( cumsum(p360, nc_dim) + a0 ); 
a0(:, 2:nc)= af(:, 1:nc-1); 

% check for zeros
test_mat= [zeros(np, 1) af];
zpos=  diff(test_mat, [], nc_dim) == 0; 
zpres= any(any(zpos));

% compute theta to evaluate arc for max non-zero percentage
theta= cell(np, nc); x= cell(np, nc);  y= cell(np, nc); inc= zeros(np, 1); 
for n= 1:np
    inz{n}= find(~zpos(n, :)); 
    [~, loc]= max(data(n, inz{n})); % index the max non-zero comp
    ii(n)= inz{n}(loc);             % to prevent too small increment
    theta{n, ii(n)}= linspace(a0(n, ii(n)), af(n, ii(n)), 300); 
    theta{n, ii(n)}= theta{n, ii(n)}(2:end);
    inc(n, 1)=   mean(diff(theta{n, ii(n)}));
end

% compute theta for remaining percentages using same increment
for n= 1:np
    rest= find(~ismember(1:nc, ii(n)));
    for p= rest
        theta{n, p}= a0(n, p)+inc(n):inc(n):af(n, p);
    end
end

if zpres    % handle case where zeros are present
    [zr, zc]= find(zpos);     % index zeros
    % if data already sums to 1 (i.e., one component accounts for 100%),
    % then complete the circle
    % else, skip the percentage, producing a partial pie
    for z= 1:length(zr)
        if sum(data(zr(z), :)) == 1    
            if zc(z) > ii(zr(z))
                theta{zr(z), zc(z)-1}(end+1)= theta{zr(z), zc(z)-1}(end)+inc(zr(z)); 
            else
                theta{zr(z), zc(z)+1}= [theta{zr(z), zc(z)+1}(1)-inc(zr(z)), theta{zr(z), zc(z)+1}];
            end
        end
    end
end

% compute circular arcs
for n= 1:np
    x(n, :)= cellfun(@(t) r*cos(t) + h(n), theta(n, :), 'UniformOutput', false); 
    y(n, :)= cellfun(@(t) r*sin(t) + k(n), theta(n, :), 'UniformOutput', false);
end

% pre-allocate structure of graphics objects for each series
for s= np:-1:1
    arcs(s).series= gobjects(1, nc);
    lbls(s).series= gobjects(1, nc); 
end

% pre-allocate text coordinates & determine position
xc=  zeros(np, nc);   yc= xc; 
centerTheta= cellfun(@median, theta); 
[xt, yt]= pol2cart(centerTheta, r*scF); 

% plot arcs & percentages
for n= 1:np
    % assign color matrix based on rule
    if strcmpi(sc, 'series')
        colors= repmat(fc(n, :), nc, 1) .* dark; 
    elseif strcmpi(sc, 'category')
        colors= fc;
    end
    % plot elements, skipping zeros (empty cells)
    for j= 1:nc
        if zpres
            if isempty(theta{n, j})
                continue
            end
        end

        xc(n, j)= xt(n, j) + h(n);
        yc(n, j)= yt(n, j) + k(n); 
        [halign, valign]= getAlignmentFromAngle(centerTheta(n, j)); % from pie.m

        arcs(n).series(1, j)= plot(x{n, j}, y{n, j}, '-', ...
                                   'Color', colors(j, :), ...
                                   'LineWidth', r*(12/(np/2)));   hold on
        lbls(n).series(1, j)= text(xc(n, j), yc(n, j), txt(n, j), ...
                                   'HorizontalAlignment', halign, ...
                                   'VerticalAlignment', valign);
    end
end
set(gcf, 'color', 'w')
axis([x_lo x_hi y_lo y_hi]) 
if np == 1
    axis square
else
    axis image; 
end
axis off; 

H.arcH=  arcs;
H.txtH=  lbls; 
H.color= colors; 

end


%% Helper functions

function [halign, valign] = getAlignmentFromAngle(angle)
% Determine the text label alignment based on the angle around the circle.

% Convert the angle to degrees
angle = (180/pi)*angle;

% Round the angles to the nearest 45 degrees
angle = mod(round(angle/45)*45,360);

% Determine the horizontal alignment
if angle == 90 || angle == 270
    halign = 'center';
elseif angle > 90 && angle < 270
    halign = 'right';
else
    halign = 'left';
end

% Determine the vertical alignment
if angle == 0 || angle == 180
    valign = 'middle';
elseif angle > 180 && angle < 360
    valign = 'top';
else
    valign = 'bottom';
end

end
