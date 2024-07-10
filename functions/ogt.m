% SUMMARY
% Created by Connor Gallimore, 11/09/2022
% 
% OWT Orange-White-Teal color map. This is a diverging colormap, 
%   useful for showing data with: 
%   1) positive values in orange, 
%   2) negative values in teal, 
%   3) and middle values in white.
%
%   OWT(M) returns an M-by-3 matrix containing a colormap. 
%   The colors begin with dark teal and increase through light teal to 
%   white, then increase to bright and then dark orange. 
%
%   OWT returns a colormap with the same number of colors
%   as the current figure's colormap. If no figure exists, MATLAB uses
%   the length of the default colormap.
%
%   EXAMPLE
%
%   This example shows how to reset the colormap of the current figure.
%
%       colormap(owt)
%
%   See also PARULA, AUTUMN, BONE, COLORCUBE, COOL, COPPER, FLAG, GRAY,
%   HOT, HSV, JET, LINES, PINK, PRISM, SPRING, SUMMER, WHITE, WINTER,
%   COLORMAP, RGBPLOT.
%--------------------------------------------------------------------------


function map = ogt(m)

if nargin < 1, m = size(get(gcf,'colormap'), 1); end


values = flipud( [0.768627450980392	 0.258823529411765	0.101960784313725;
                  0.976470588235294	 0.560784313725490	0.270588235294118;
                  0.996078431372549	 0.749019607843137	0.431372549019608;
                          0.9              0.9                 0.9;
                  0.592156862745098	 0.807843137254902	0.800000000000000;
                  0.070588235294117  0.564705882352941	0.556862745098039;
                  0.086274509803922	 0.349019607843137	0.290196078431373] );

P = size(values, 1);
map = interp1(1:P, values, linspace(1, P, m), 'linear');
