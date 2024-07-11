%SUMMARY
% Author: Connor Gallimore
% 12/21/2022

% Creates clock labels using either regular (AM/PM) or military time on a
% 24-hour day.

% Required arguments:
    % 't_min', set as 0 for a range beginning at midnight
    % 't_max', set as 24 to consider a full 24-hour day, back to midnight
    % 'n_hrs', a positive integer denoting the desired number of HH:mm ticks
    % 'time_type', set as character vectors 'regular' or 'military'

% Example:
% labels = makeClockLabels(0, 24, 12, 'regular')
% outputs a 1 x n_hrs cell of character vector labels in 2-hour increments
% for a full 24-hour day, including AM/PMs appended to applicable hour
% marks

% Note: if the requested time_type is 'regular', AND the requested number 
% of ticks results in exact hour times, the function chops the minutes and 
% simply returns the hour + AM/PM
%--------------------------------------------------------------------------


function labels = makeClockLabels(t_min, t_max, n_hrs, time_type)

expectedTimeTypes= {'military', 'universal', '24H', '24Hr', ...
                    'civilian', 'regular', '12H', '12Hr'};

tt_input= find(cellfun(@(c) strcmpi(c, time_type), expectedTimeTypes));

label_tmp= t_min: (t_max / n_hrs) : t_max;
label_tmp= label_tmp(1:end-1); 
min_str_fin= repmat("00", 1, length(label_tmp)); 
clock_labls= strings(1, length(min_str_fin)); 
hr= floor(label_tmp); 

if any(mod(label_tmp, 1) ~= 0)
    loc_NI= mod(label_tmp, 1) ~= 0;            % idx non-integers
    nonInts= label_tmp(loc_NI);                % grab non-integers
    min_of_hr= round(mod(nonInts, 1) * 60);    % convert to min of hour

    % if any of the minute digits are single, add leading zero
    min_digits= ceil(log10(max(1, abs(min_of_hr) * (1+eps))));
    min_digits(min_digits == 0) = 1; 
    one_dig= min_digits == 1;
    min_str= string(min_of_hr); 
    min_str(one_dig) = strcat("0", min_str(one_dig)); 
    min_str_fin(loc_NI)= min_str; 
else
    min_of_hr= zeros(1, length(label_tmp));
end

if ismember(tt_input, 1:4)  % strcmpi(time_type, 'military')
    % if any of the hour digits are single, add leading zero
    hr_digits= ceil(log10(max(1, abs(hr) + 1)));  
    hr_digits(hr_digits == 0) = 1; 
    one_dig= hr_digits == 1;
    clock_labls(one_dig) = strcat("0", string(hr(one_dig)), ":", min_str_fin(one_dig)); 
    clock_labls(~one_dig) = strcat(string(hr(~one_dig)), ":", min_str_fin(~one_dig)); 

else
    % add AMs and PMs, dont add leading zeros for regular time
    AMPM= min_str_fin; 
    afternoon=  label_tmp >= 12;
    afternoon2= label_tmp > 12;
    AMPM(~afternoon) = "AM"; AMPM(afternoon)= "PM"; 
    reg_labls= [hr(~afternoon2), hr(afternoon2) - 12];
    reg_labls(1)= 12; 

    % if increments result in all exact hour marks, chop off minutes
    if any(min_of_hr)
        clock_labls= strcat(string(reg_labls), ":", min_str_fin, AMPM); 
    elseif all(~min_of_hr)
        clock_labls= strcat(string(reg_labls), AMPM); 
    end

end

labels= cellstr(clock_labls); 

end