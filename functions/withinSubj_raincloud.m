%SUMMARY
% Connor Gallimore
% 01/31/2023

% Within-subjects raincloud, requires n x 2 column matrix of data where
% each condition has an equal number (n) of observations

% Modified from: rm_raincloud
% https://github.com/RainCloudPlots/RainCloudPlots/blob/master/tutorial_matlab/rm_raincloud.m
%--------------------------------------------------------------------------

function h = withinSubj_raincloud(data, colors, varargin)

% process optional Name,Val pairs, assign defaults if absent
tmpO= strcmpi(varargin, 'orientation'); 
tmpB= strcmpi(varargin, 'bandwidth'); 
tmpA= strcmpi(varargin, 'alpha'); 
tmpD= strcmpi(varargin, 'density');
tmpR= strcmpi(varargin, 'dropSize');
tmpL= strcmpi(varargin, 'lineColor');

if any(tmpO); plot_ori= varargin{find(tmpO) + 1}; 
else;         plot_ori= 0;      % left-to-right plotting by default
end

if any(tmpB); bandwidth= varargin{find(tmpB) + 1}; 
else;         bandwidth= [];    % let the function specify the bandwidth
end

if any(tmpA); alpha= varargin{find(tmpA) + 1}; 
else;         alpha= 0.35;      % default transparency of plot elements
end

if any(tmpD); density= varargin{find(tmpD) + 1}; 
else;         density= 200;     % default # of probability density bins
end

if any(tmpR); raindrop_Sz= varargin{find(tmpR) + 1}; 
else;         raindrop_Sz= 24;  % default raindrop size, NOTE: was 60
end

if any(tmpL); lineColor= varargin{find(tmpL) + 1}; 
else;         lineColor= [0.5 0.5 0.5];  % default line colors
end

% pre-allocate properties of density plots
[~, n_meas] = size(data);
n_bins = repmat(density, 1, n_meas);
ks= cell(1, n_meas); x= ks; 

for n= 1:n_meas
    % compute density using 'ksdensity'
    [k, y, ~] = ksdensity(data(:, n), 'NumPoints', n_bins(n), 'bandwidth', bandwidth);
    ks{n}=k;  x{n}=y;

    q{n}     = (1:n_bins(n) - 1)';
    faces{n} = [q{n}, q{n} + 1, q{n} + n_bins(n) + 1, q{n} + n_bins(n)];
end

spacing     = 2 * mean(mean(cellfun(@max, ks)));
ks_offsets  = (0:n_meas-1) .* spacing;

% flip so first plot in series is plotted on the *top*
ks_offsets  = fliplr(ks_offsets);

for n = 1:n_meas
    verts{n} = [x{n}', ks{n}' + ks_offsets(n); x{n}', ones(n_bins(n), 1) * ks_offsets(n)];
end

% jitter for the raindrops
jit_width = spacing / 6;   % default 8

for n = 1:n_meas
    jit{n} = jit_width + rand(length(data(:, n)), 1) * jit_width;
end

summaryAVS=  mean(data, 1, 'omitnan'); 

hold on

% plot stuff in series
for n = 1:n_meas    
    % plot patches
    h.cloud_patch(n) = patch('Faces', faces{n}, 'Vertices', verts{n}, 'FaceVertexCData', colors(n, :), ...
                             'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', alpha);
    
    if n == n_meas
        patch_Y = get(h.cloud_patch(n), 'YData');
        new_Y =  -1 * patch_Y' ;
        set(h.cloud_patch(n), 'YData', new_Y);
        jit{n}= -jit{n};
    end   

    % scatter raindrops
    h.raindrops(n) = scatter(data(:, n), -jit{n} + ks_offsets(n), 'MarkerFaceColor', colors(n, :), ...
                             'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha, 'SizeData', raindrop_Sz);

    if n == n_meas
        xC= [get(h.raindrops(n-1), 'XData'); get(h.raindrops(n), 'XData')];
        yC= [get(h.raindrops(n-1), 'YData'); get(h.raindrops(n), 'YData')];
        h.ws_lines(:, n-1)= plot(xC, yC, '-', 'color', [lineColor, alpha]);  
    end
end

% mean line
h.mean_line(1) = line(summaryAVS(1:2), ks_offsets(1:2), 'LineWidth', 4, 'Color', lineColor);

% mean dots
for n = 1:n_meas    
    h.mean_dots(n) = scatter(summaryAVS(n), ks_offsets(n), 'MarkerFaceColor', colors(n, :), ...
                             'MarkerEdgeColor', [0 0 0], 'SizeData', raindrop_Sz * 2, 'LineWidth', 2);
end

set(gca, 'YTick', fliplr(ks_offsets), 'YTickLabel', n_meas:-1:1);

% determine plot orientation
% default is left-to-right, plot_ori can be set to 1 
% NOTE: Because it's easier, we actually prepare everything plotted
% top-to-bottom, then - by default - we rotate it here. 

if ~plot_ori
    view([90 -90]);
    axis ij
end

end