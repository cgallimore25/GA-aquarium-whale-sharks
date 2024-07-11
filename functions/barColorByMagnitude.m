%SUMMARY
% Author: Connor Gallimore
% 01/25/2023
%--------------------------------------------------------------------------


function bh = barColorByMagnitude(y_data, colorMap, varargin)


edges= linspace(min(y_data), max(y_data), length(colorMap)); 
bin_obs= discretize(y_data, edges); 

sz= length(y_data); 

barH=  gobjects(1, sz);
xtips= zeros(1, sz); 
ytips= zeros(1, sz);

for k= 1:length(y_data)
    barH(k)= barh(k, y_data(k));  hold on; 
    set(barH(k), 'FaceColor', colorMap(bin_obs(k), :))
    xtips(k)= barH(k).YEndPoints + 0.3;  
    ytips(k)= barH(k).XEndPoints; 
end

bh.barH= barH; 
bh.xtips= xtips;
bh.ytips= ytips; 

end