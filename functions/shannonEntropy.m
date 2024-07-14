function H = shannonEntropy(p, varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% matrix assumes a row matrix where rows are different collections of
% probabilities and columns represent the individual probabilities in each
% collection

% zeros will be returned as NaNs

% identify input type
matrix= @(x) isnumeric(x) && ~isscalar(x) && ndims(x) >= 2 && ~isvector(x);
vector= @(x) isnumeric(x) && ~isscalar(x) && isvector(x);

dim_arg= @(x) isnumeric(x) && isscalar(x) && x > 0;

dim_tmp= cellfun(dim_arg, varargin); 

if vector(p)
    H= -sum(p .* log2(p), "omitnan"); 
elseif matrix(p)
    if any(dim_tmp)
        dim= cell2mat(varargin);
    else
        dim= 2; 
    end
    H= -sum(p .* log2(p), dim, "omitnan"); 
end

end