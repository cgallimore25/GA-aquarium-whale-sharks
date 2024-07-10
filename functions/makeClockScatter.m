%SUMMARY
% Author: Connor Gallimore
% 12/21/2022

% This function creates a polar scatter in the form of a clock plotting a
% datetime array, useful for visualizing circadian data. against some 
% distance variable. 
% By default, the function assumes military time, unless specified 
% otherwise (see Optional arguments). 

% Required arguments:
    % 'datetime_data', a datetime array type (see MATLAB documentation)
    % 'rho', some distance variable where higher values lie further away
    %        from the center
    % 't_min', time to set at top of clock, usually 0
    % 't_max', end time that would "close" interval (overlap with 0), 
    %          e.g. 24 for a clock plotting a 24-hr day
    % 'n_hrs', a positive integer denoting the desired number of HH:mm ticks

% Optional arguments:
% should be entered in the form of Name,Value pairs
    % 'DotSize', a positive integer indicating size of scatter dots
    % 'Color', supported line/marker spec or RGB triplet
    % 'EdgeColor', supported line/marker spec or RGB triplet
    % 'timeMethod', character vector specifying 'regular' (AM/PM) or 
    %               'military'
    % 'inclHand', adds an arrow at an additionally specified single point
    %             across the full distance of rho

% Outputs
% 1 output argument returns a structure of scatter and PolarAxes handles
% 2 output args adds the non-normalized times in the form of whole hour and
%               fraction the hour following the decimal (double)
% 3 output args adds the times normalized to the 2pi polar axes
% 4 output args returns the clock labels created by makeClockLabels(), see
%               its independent function file for more on this
%--------------------------------------------------------------------------


function varargout = makeClockScatter(datetime_data, rho, t_min, t_max, n_hrs, varargin)

% search for possible Name,Val pairs 
size_tmp=  strcmpi(varargin, 'DotSize'); 
face_col=  strcmpi(varargin, 'Color'); 
edge_col=  strcmpi(varargin, 'EdgeColor'); 
tMeth_tmp= strcmpi(varargin, 'timeMethod'); 
hand_tmp=  strcmpi(varargin, 'inclHand');

% process optional Name,Val pairs and set defaults if none specified
if any(size_tmp)
    dotSz= varargin{find(size_tmp) + 1};
else
    dotSz= 48;
end

if any(edge_col)
    ec= varargin{find(edge_col) + 1};
else
    ec= 'k';
end

if any(tMeth_tmp)
    time_type= varargin{find(tMeth_tmp) + 1};
else
    time_type= 'military';
end

degs= linspace(0, 360-(360/n_hrs), n_hrs); 

% separate fxn converting hours since midnight to units of seconds,
% and normalizing to clock interval
[t_norm, clock_times]= normalizeTimeCircadianClock(datetime_data, t_min, t_max);

% separate fxn that makes equally spaced clock labels
labels= makeClockLabels(t_min, t_max, n_hrs, time_type);

if any(face_col)
    fc= varargin{find(face_col) + 1};
    if strcmpi(fc, 'byDensity')
        [~, t_norm, rho, fc]= scatter_kde(t_norm, rho, 'DenseSort', 'ValsOnly');
    end
else
    fc= 'c';
end

% make plot and axis handles
h= polarscatter(t_norm, rho, dotSz, fc, 'filled', 'MarkerEdgeColor', ec); hold on; 
ax= gca;

% plot hand at theta coordinate, if any specified by 'incMedHand'
if any(hand_tmp)
    H_pos= varargin{find(hand_tmp) + 1};
    [norm_th, ~]= normalizeTimeCircadianClock(H_pos, t_min, t_max);
    rsz= size(rho);  mh_theta= repmat(norm_th, rsz); 
    clockHand= polarplot(mh_theta, rho, '-k'); 
end

% orient time clockwise with midnight at the top
ax.ThetaDir= 'clockwise';  
ax.ThetaZeroLocation = 'top';

% change labels to custom labels
thetaticks(degs); 
thetaticklabels(labels); 

% pack plot elements into struct
clockscatter.handle= h;
clockscatter.axes= ax;
clockscatter.thetalabels= ax;

if exist("clockHand", "var")
    clockscatter.Hand= clockHand;
end

if any(face_col)
    if strcmpi(varargin{find(face_col) + 1}, 'byDensity')
        clockscatter.colors= fc; 
    end
end

% outputs
varargout= { clockscatter, clock_times, t_norm, labels };

end