%SUMMARY
% Author: Connor Gallimore
% 12/21/2022

% This function creates a polar histogram in the form of a clock plotting a
% datetime array, useful for visualizing circadian data. By default, the 
% function assumes military time, unless specified otherwise 
% (see Optional arguments). 

% Required arguments:
    % 'datetime_data', a datetime array type (see MATLAB documentation)

% Optional positional arguments:
    % t_min - time to place at top of clock (usually 0 for midnight); can
    %         input [] to use default and access later positionals 
    % t_max - max time on interval, 'closes' polar data back to 0; can
    %         input [] to use default and access third positional n_hrs
    % n_hrs - positive scalar integer, number of tick increments

    % Examples:
    % makeClockHist(datetime_data, t_min, t_max, n_hrs, nbins) uses the
    % number of bins specified by the positive integer, nbins

    % makeClockHist(datetime_data, t_min, t_max, n_hrs, edges) sorts
    % datetime_data into bins with bin edges specified by the vector, edges

    % makeClockHist(datetime_data, [], [], n_hrs, edges) defaults first 2
    % and places datetime_data into bins with bin edges specified by the 
    % vector, edges, displaying specified number of hours as ticks

% Name,Value pairs
    % 'FaceColor', supported line/marker spec or RGB triplet
    % 'FaceAlpha', scalar from 0-1, 1 = completely opaque, 0 = transparent
    % 'EdgeColor', supported line/marker spec or RGB triplet
    % 'timeMethod', character vector specifying 'regular' (AM/PM) or 
    %               'military'
    % 'Normalization', supported 'count' or 'probability' inputs

% Outputs
% 1 output argument returns a structure of histogram and PolarAxes handles
% 2 output args adds the non-normalized times in the form of whole hour and
%               fraction the hour following the decimal (double)
% 3 output args adds the times normalized to the 2pi polar axes
% 4 output args returns the clock labels created by makeClockLabels(), see
%               its independent function file for more on this
%--------------------------------------------------------------------------


function varargout = makeClockHist(datetime_data, varargin)

% set defaults
defaultT0= 0;           % time at top
defaultTN= 24;          % close data with end time
defaultHr= 12;          % hour increments
defaultnB= 51;          % number of bins
defaultFC= 'c';         % histogram face color
defaultFA= 0.6;         % histogram face alpha
defaultEC= 'k';         % histogram edge color
defaultTT= 'military';  % time-type
defaultN=  'count';     % normalization method

expectedTimeTypes= {'military', 'universal', '24H', '24Hr', ...
                    'civilian', 'regular', '12H', '12Hr'};
expectedNormMethod= {'count', 'probability'};

p= inputParser;         % inputParser object

% make some validation functions
scalarNum= @(x) isnumeric(x) && isscalar(x);   % to avoid repetition
validRGB=  @(x) all(size(x) == [1 3]) && all(x <= 1) && all(x >= 0);
validData= @(x) isdatetime(x);
validT0_TN= @(x) scalarNum(x) || isempty(x);
validNumHrs= @(x) scalarNum(x) && (x >= 1) && ~mod(x, 1);
checkBinNum= @(x) scalarNum(x) && (x > 1);
checkBinEdg= @(x) ~isscalar(x) && isvector(x) && isnumeric(x);
validBinMethod= @(x) checkBinNum(x) || checkBinEdg(x);
validColor= @(x) ischar(x) || validRGB(x);
validAlpha= @(x) scalarNum(x) && (x >= 0) && (x <= 1);

% add all inputs
addRequired(p, 'datetime_data', validData); 
addOptional(p, 't_min', defaultT0, validT0_TN);
addOptional(p, 't_max', defaultTN, validT0_TN); 
addOptional(p, 'n_hrs', defaultHr, validNumHrs); 
addOptional(p, 'binMethod', defaultnB, validBinMethod); 
addParameter(p, 'FaceColor', defaultFC, validColor); 
addParameter(p, 'FaceAlpha', defaultFA, validAlpha);
addParameter(p, 'EdgeColor', defaultEC, validColor); 
addParameter(p, 'timeMethod', defaultTT, @(x) any(validatestring(x, expectedTimeTypes))); 
addParameter(p, 'Normalization', defaultN, @(x) any(validatestring(x, expectedNormMethod))); 

% parse and assign all inputs
parse(p, datetime_data, varargin{:});

fc= p.Results.FaceColor;
fa= p.Results.FaceAlpha;
ec= p.Results.EdgeColor;
t_min= p.Results.t_min;
t_max= p.Results.t_max; 
n_hrs= p.Results.n_hrs;
time_type= p.Results.timeMethod;
h_tmp= p.Results.binMethod; 
norm_mth= p.Results.Normalization; 

if isempty(t_min); t_min= defaultT0; end  % if empty, replace w/ defaults
if isempty(t_max); t_max= defaultTN; end

if checkBinEdg(h_tmp)  % if vector, normalize to polar coords
    hist_arg= (h_tmp - min(h_tmp)) / (max(h_tmp) - min(h_tmp)) * 2 * pi; 
else
    hist_arg= h_tmp;   % else, must be bin number
end

degs= linspace(0, 360-(360/n_hrs), n_hrs); % ticks plotted (later replaced)

% convert hours since midnight to units of seconds, normalize to clock interval
[t_norm, clock_times]= normalizeTime2Polar(datetime_data, t_min, t_max);

% separate fxn that makes equally spaced clock labels
labels= makeClockLabels(t_min, t_max, n_hrs, time_type);

% make plot and axis handles
h= polarhistogram(t_norm, hist_arg, 'FaceColor', fc, 'EdgeColor', ec, 'FaceAlpha', fa, 'Normalization', norm_mth); 
ax= gca;

% orient time clockwise with midnight at the top
ax.ThetaDir= 'clockwise';  
ax.ThetaZeroLocation = 'top';

% change labels to custom labels
thetaticks(degs); 
thetaticklabels(labels); 

% pack plot elements into struct
clockhist.handle= h;
clockhist.axes= ax;
clockhist.thetalabels= labels;
clockhist.bin_arg= hist_arg; 

% outputs
varargout= { clockhist, clock_times, t_norm, labels };

end


%% Local funcs

function varargout = normalizeTime2Polar(datetime_data, t_min, t_max)

% convert hours since midnight to units of seconds
secs_since_midnight= seconds( timeofday(datetime_data) ) / 3600;

% normalize it to the clock interval
t_norm= (secs_since_midnight - t_min) / (t_max - t_min) * 2 * pi; 

varargout= {t_norm, secs_since_midnight}; 

end

