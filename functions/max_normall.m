% SUMMARY
%  Author: Connor Gallimore

%  This function re-scales all values of a matrix from 0-1, such that the
%  lowest value in the matrix (e.g. min(m(:)) ) becomes 0, and the highest 
%  value becomes 1. 
%--------------------------------------------------------------------------


function [y] = max_normall(x)

y = x - min(x(:));
y = y ./ max(y(:));
    
end

