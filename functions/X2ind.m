function s = X2ind(Obs_mat, Ex_mat, alpha)
% Chi-squared GOF test

n_rows= size(Obs_mat, 1);
n_cols= size(Obs_mat, 2);

Res= Obs_mat - Ex_mat;   % Residual -- difference between predicted values and observed values
SSres= Res.^2;           % Residual sum of squares -- amount of variance unexplained by model
comp=  SSres ./ Ex_mat;
X2= sum( comp );

df= (n_rows-1) * (n_cols-1);       % degrees of freedom -- number of groups minus 1

X2crit= chi2inv(1-alpha, df);

s.Resid= Res;  
s.SSres= SSres;
s.Computed= comp;
s.X2= X2;  
s.df= df;
s.X2crit= X2crit;

end

