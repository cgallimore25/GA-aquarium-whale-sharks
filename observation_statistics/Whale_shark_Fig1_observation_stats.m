%% Whale shark observation statistics


%% Parameters

% make some colormaps
cmap= cmocean('matter');  cm= flipud(cmap); 

m_owt= [0.768627450980392	 0.258823529411765	0.101960784313725;
        0.976470588235294	 0.560784313725490	0.270588235294118;
        0.996078431372549	 0.749019607843137	0.431372549019608;
                1               1                   1;
        0.592156862745098	 0.807843137254902	0.800000000000000;
        0.070588235294117    0.564705882352941	0.556862745098039;
        0.086274509803922	 0.349019607843137	0.290196078431373]; 
    
tealDark= m_owt(6, :);  tealLight= m_owt(5, :); 


% make this your path name where all whale shark stuff is kept
saveDir= pwd;


%% Load data

fileName= 'GAA_whale_shark_master_datasheet.csv'; 

opts= detectImportOptions( fileName );

opts= setvartype(opts, {'Day', 'Observer', 'Location', ...
                        'Shark', 'Depth', 'Direction_of_Swim', ...
                        'Date', 'Time_of_Day', 'Away_from_Wall', ...
                        'Other_Notes'}, 'string');      

dataTbl = readtable(fileName, opts);   


%% Process datetime data

dates= dataTbl.Date;        missing_dates=   find(ismissing(dates)); 
dates_dt= datetime(dates, 'InputFormat', 'MM/dd/uuuu', 'TimeZone','America/New_York');
TOD= dataTbl.Time_of_Day;   missing_tstamps= find(ismissing(TOD)); 
TOD_dt=   datetime(TOD, 'InputFormat', 'HH:mm');
day_col=  categorical( lower( dataTbl.Day ) );  
month_col= dataTbl.Month;  [~, mo_ix]= sort(unique(month_col, 'stable')); 

datetimes_all= dates_dt + timeofday(TOD_dt);  


%% Process month names and make labels
month_names= month(datetimes_all, 'name');  ix_tmp= find(~cellfun(@isempty, month_names)); 
month_str= string( cellfun(@(s) s(1:3), month_names(ix_tmp), 'UniformOutput', false) );

month_labels_long= unique( string(month_names(ix_tmp)) , 'stable');  month_labels_long= month_labels_long(mo_ix); 
month_labels_short= unique(month_str, 'stable');  month_labels_short= month_labels_short(mo_ix); 


%% Plot Figure 1 Panels C-F -- Observation statistics

% Time histograms by all, month, year, and day show that overwhelming 
% majority of data is taken from Mar-May, most notes were taken Mon-Thurs,
% and 2009 was the year with most notes

yrs_active= unique(datetimes_all.Year); 
yrs_active= yrs_active(~isnan(yrs_active)); 
ya= datetimes_all.Year; 

n_obs= length(datetimes_all); 

C_title= strcat('All timestamped observations', {' '}, '(N=', num2str(n_obs), ')'); 

figure('Name', 'Timestamp histograms'); 
subplot(6, 2, 1:2:3);    histogram(datetimes_all, 'FaceColor', tealDark, 'FaceAlpha', 1);  
                         ylabel('Count'); grid on; 
                         title(C_title)  % set(gca, 'GridAlpha', 1); 
subplot(6, 2, 2:2:4);    th= stackedHist(month_col, findgroups(ya), string(yrs_active));
                         grid on; colormap(cmap);  % set(gca, 'GridAlpha', 1); 
                         xticks(1:12);  xticklabels(month_labels_short); xlim([0 13]); title('Month')
subplot(6, 2, 7:2:11);   [tph, clock_times]= makeClockHist(datetimes_all, 0, 24, 12, 'FaceColor', tealLight, 'EdgeColor', 'none', 'FaceAlpha', 1); 
                         title('Time of day');  % set(gca, 'GridAlpha', 1);  % default gridalpha 0.15

ct_vec= histcounts(clock_times, 0:1:24);
pct_t_obs= ct_vec./sum(ct_vec) * 100;

% get some day-of-week statistics to finish plot
day_ord= reordercats(day_col, [2 6 7 5 1 3 4]); 
dn= string( cellfun(@(s) s(1:3), categories(day_ord), 'UniformOutput', false) );
day_cnts= histcounts(day_ord);
mon_cnts= histcounts(month_col);  
mon_pcts= mon_cnts ./ sum(mon_cnts) * 100;                         
                         
subplot(6, 2, 8:2:12);   bd= barColorByMagnitude(fliplr(day_cnts), cmap); grid on; % set(gca, 'GridAlpha', 1); 
                         yticks(1:length(dn));  
                         yticklabels( flipud(dn) );
                         set(gca, 'xlim', [0 8000], 'XAxisLocation', 'top'); 
                         set(gcf, 'color', 'w');
                         text(bd.xtips' + 100, bd.ytips', string(fliplr(day_cnts))', 'VerticalAlignment', 'middle')
                         title('By weekday')


%% Plot Figure 1 Panel G -- Location of observer

Location = categorical( dataTbl.Location ); 

Loc_cat= reordercats(Location, [2 3 1 6 5 7 8 4]);
[loc_g, L_cats]= grp2idx(Loc_cat); 
n_locs= length(L_cats); 
inc=    360/n_locs; 
L_num= 0:8;  
rhoL= 3 * ones(1, length(L_cats)); 

L_cnts= histcounts(Loc_cat); 
L_pct= L_cnts / sum(L_cnts) * 100; 
L_norm= (L_num - min(L_num)) / (max(L_num) - min(L_num)) * 2 * pi; 

% observation location table
obs_loc_tbl= table(L_cats, L_cnts', L_pct', 'VariableNames', {'Location', 'Counts', '% Total'});

figure('Name', 'Ocean Voyager observer location'); 
loc_pb= polarbubblechart(L_norm(1:8), rhoL, round(L_pct, 2), 'r', 'MarkerFaceAlpha', 1);
b_sz= [10 50];
bubblesize(b_sz); 
bubblelim([min(round(L_pct, 2)) max(round(L_pct, 2))])
bubblelegend('Observation density', 'Location', 'eastoutside')
thetaticks(0:inc:(360-inc))
thetaticklabels(L_cats)
set(gca, 'ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top'); 

