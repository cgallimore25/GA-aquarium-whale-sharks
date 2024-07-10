% SUMMARY
% Author: Connor Gallimore
% 01/19/2023

% Normalize datetime data to a user-defined min (e.g. 0) and max (e.g. 24) 
% hours on a clock interval
%--------------------------------------------------------------------------


function varargout = normalizeTimeCircadianClock(datetime_data, t_min, t_max)

% convert hours since midnight to units of seconds
secs_since_midnight= seconds( timeofday(datetime_data) ) / 3600;

% normalize it to the clock interval
t_norm= (secs_since_midnight - t_min) / (t_max - t_min) * 2 * pi; 

varargout= {t_norm, secs_since_midnight}; 

end