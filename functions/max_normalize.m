%SUMMARY
% Author: Connor Gallimore
% 10/26/2022

% Normalize a N x M matrix by scaling each row (N) from 0 (min value) to 1 (max value)

% Written for use on time-series data where each row contains the nth
%   series and each column contains the mth time point 
%   (e.g. neuronal recordings of neurons (rows) by time (cols) )
%--------------------------------------------------------------------------

function y = max_normalize(x)

row_min = min(x, [], 2);
y = x - row_min;
row_max = max(y, [], 2);
y= y ./ row_max;

% for loop version (slower)
% y= zeros( size(x) );
%     for r= 1:size(x, 1)
%         y(r, :)= x(r, :) - min( x(r, :) );
%         y(r, :)= y(r, :) ./ max( y(r, :) );
%     end

end

