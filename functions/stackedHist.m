% Author: Connor Gallimore
% 04/02/2022

% SUMMARY
%   Plot a stacked histogram where the order of stacking, and thus the color gradient, 
%   can be specified by the user (see Op Name,Val pairs) based on some metric computed 
%   for each group.  The default behavior assumes the distribution is discrete 
%   (however, see optional argument 'binRange' to treat your data as continuous). 

%   One of my motivations for writing this was that the default behavior of 
%   bar(data, 'stacked') shows the bottom bar and corresponding color at the 
%   top of the legend, such that the order is visually flipped. This kind
%   of irked me. While less efficient, this code orders them correctly by
%   manually building the stacked bars from highest magnitude to lowest,
%   and assigns CData for each category so that a colormap() call may follow
%   in your main script.

% Required:
%   data        -- 1 x n vector of all data points to plot in the histogram
%   groupidVec  -- 1 x n vector (equal in size to 'data') with group identity
%                  coded as an arbitrary number for each group
%   groups      -- 1 x n string array of group names (used to deduce the
%                  number of groups

% Usage
%   stackedHist(data, groupidVec, groups)

%   alternatively, if dont care about the data matrix being plotted or sorted index (spec),
%   use 'h = stackedHist(data, groupidVec, groups)' 

% Optional args
%   a positive scalar from 0-1 indicating the width of the bar, where 1 shows
%                   all bars touching (default), mimicking a histogram
%   'probability' - this tells the function you want to normalize the individual categories 
%                   (rows) by the sum of their elements
%   'valsOnly' - if you only want the data matrix, varargout will truncate the function 
%                prior to plotting, returning up to 2 outputs: 'counts' (matrix), 
%                and/or corresponding 'x_vec' for the bins

% Optional Name,Value pairs (varargin)
%   'sortColorsBy' - the next input you pass should be a 1 x n vector (one element per group) 
%                    of some metric to stack the histogram and color by. You can then
%                    follow with a sorting method 'ascend' or 'descend'
%   'binRange' - if you have a continous distribution (>50 points) and/or you want to bin over
%                a range of values that DO NOT begin at 0 and/or DO NOT increment by 1, you 
%                can pass a custom bin range ( e.g. 6.0 : 0.1 : 8.0; )
%   'EdgeColor' - controls color of surrounding line, typical character args, or 'none', apply

% Optional outputs (varargout)
%   1 output argument will supply a structure with a handle to each bar, the legend, 
%     and a structure to color data tick marks
%   2 output arguments will supply the above in addition to the matrix being plotted
%   3 output arguments will supply the above in addition to the order of plotting, 
%     if you chose to sort them in anyway with optional input 'sortColorsBy'
%--------------------------------------------------------------------------

function varargout = stackedHist(data, groupidVec, groups, varargin)

    % Define & pre-allocate some vars
    nG=  length(groups);    lgdStrs=  cell(1, nG);    

    % Parse optional inputs
    tmpS= strcmpi(varargin, 'sortColorsBy');        % sort order
    tmpB= strcmpi(varargin, 'binRange');            % bin range
    tmpN= strcmpi(varargin, 'probability');         % normalization
    tmpE= strcmpi(varargin, 'EdgeColor');           % edge arg (bar)
    tmpW= isscalar(varargin) & isnumeric(varargin); % width scalar
    tmpV= strcmpi(varargin, 'valsOnly');            % if dont want plot
    
    % Assign width (bars)
    if any(tmpW)
        width= varargin{tmpW};
    else
        width= 1; % default width of 1 (i.e. touching)
    end
    
    % Assign sort order
    if any(tmpS)
        sOrder = varargin{find(tmpS) + 1};
        if any( [strcmpi(varargin, 'ascend'), strcmpi(varargin, 'descend')] )
            sMethod = varargin{find(tmpS) + 2};
            [~, ix]= sort(sOrder, sMethod);
        else
            [~, ix]= sort(sOrder, 'ascend') ;                    % default 'ascend' method        
        end
    else 
        ix= 1:nG;                                                % if no method specified, plot groups as arranged
    end

    group_order= ix;
    
    % Assign bin range
    if any(tmpB)
        bin_rng = varargin{find(tmpB) + 1};
        samp_rate= diff(bin_rng);   shift= samp_rate(1) / 2;     % if bin range specified, place bar in center of bin edges
    else 
        shift= 1;
        bin_rng=  0:max(data) + shift;                           % if no bin range specified, just sample from 0 : 1 : max
    end
    
    nbins= length(bin_rng)-1; 
    x_vec= bin_rng(2:end)-shift;
    counts= zeros(nG, nbins);   dat2plt= counts;

    % Create bin count vector for each Group, stack into matrix
    for g= 1:nG
        counts(g, :)= histcounts(data(groupidVec == g), bin_rng);      
    end

    if any(tmpN)   % returns probability rather than counts, if specified
        n_group_obs= sum(counts, 2);
        counts= counts ./ n_group_obs; 
    end

    if any(tmpV)   % if you only want to data matrix
        varargout= { counts, x_vec }; 
        return
    end
    
    % Build stacked hist matrix
    for g= 1:nG
        if g == 1
            dat2plt(g, :)=  counts( ix(g), :); 
        else
            dat2plt(g, :)=  dat2plt(g-1, :) + counts( ix(g), :);
        end
    end

    if any(tmpE)   % specify bar edge color
        ec= varargin{find(tmpE) + 1};
    else
        ec= 'k';   % default to black
    end

    % Plot data and output structure of handles to plot elements
    
    % Stacked bar() method that works, but exhibits the legend behavior i dislike
	% h.bars=  bar(x_vec, counts', width, 'stacked', 'FaceColor', 'flat', 'EdgeColor', ec);
    
    for b= 2:nG+1
        h.bars(b-1)= bar(x_vec, dat2plt( (-1*b) + (g+2), :), width, ...
                         'FaceColor', 'flat', 'EdgeColor', ec);  hold on;  
        h.bars(b-1).CData = (-1*b) + (g+2);  
        lgdStrs(b-1)= cellstr( groups( ix( (-1*b) + (g+2) ) ) );
    end
    clim([1 max([h.bars.CData])]);      % use full range of cmap
    set(gcf, 'color', 'w');              % tidy up
    h.lgd= legend(lgdStrs); 
    
    varargout = { h, dat2plt, x_vec, group_order };

end