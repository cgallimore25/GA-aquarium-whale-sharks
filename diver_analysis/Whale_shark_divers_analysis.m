%% Whale shark divers analysis

% Reproduces panels a-c for Figure S7 of [manuscript] 


%% Parameters

exportFig=    0;        % sets alpha to 1 if exporting raincloud plots
computeSpeed= 1;        % make 0 to pull pre-comp'd speeds from csv file
speed_units= 'meters';  % change to feet for ft/sec


if strcmpi(speed_units, 'meters')
    spd_ax_lbl= 'Speed (m/sec)'; 
    convF= unitsratio('meters', 'feet');
    spd_ci_lim= [-0.15 0.15];
else
    spd_ax_lbl= 'Speed (ft/sec)'; 
    convF= 1;   % keep units same if want in ft/sec
    spd_ci_lim= [-0.45 0.45]; 
end

% make some colormaps
m_owt= flipud( [0.768627450980392	 0.258823529411765	0.101960784313725;
                  0.976470588235294	 0.560784313725490	0.270588235294118;
                  0.996078431372549	 0.749019607843137	0.431372549019608;
                         0.9                0.9                0.9;
                  0.592156862745098	 0.807843137254902	0.800000000000000;
                  0.070588235294117  0.564705882352941	0.556862745098039;
                  0.086274509803922	 0.349019607843137	0.290196078431373] );
    
tealDark= m_owt(2, :);  tealLight= m_owt(3, :); orangeMed= m_owt(6, :); 

colors= [orangeMed; tealDark]; 

rc_col= [[164 63 97] ./ 255; [.3 .3 .3]]; 


%% Load data

fileName= '18_19_Dec_Master_data_sheet_proc2.csv'; 

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


%% Process shark data and plot some shark descriptives

Shark= upper( dataTbl.Shark );  Shark= Shark(~ismissing(Shark)); % get rid of missing shark

[Shark_g, Shark_Name]= findgroups(Shark);
n_Sharks= length(Shark_Name); 


datetimes_sharks= datetimes_all(~ismissing(Shark));  % get all datetimes where shark was noted

Dist=  dataTbl.Distance_ft;        Dist= Dist(~ismissing(Shark)); 
Time=  dataTbl.Time_sec;           Time= Time(~ismissing(Shark));

% Can either calculate speed (method 1), or pull computed speed (method 2)
if computeSpeed == 1
    Speed= (Dist * convF) ./ Time;   
    Speed(Speed > (10 * convF)) = NaN;  % nan out unreliable vals
else
    Speed= dataTbl.Speed_ftPerSec * convF;  
    Speed= Speed(~ismissing(Shark)); % get rid of missing vals
end


%% Interactions with divers

Divers= dataTbl.Divers;         Divers= Divers(~ismissing(Shark));


%% Process datetimes with / without divers present
% the below routine compares average speed times across all sharks when
% divers were or were not present in a *within-day* fashion. 
% In other words, it pulls and averages over the speed measurments for a 
% given day when divers were present, then does so when divers were present 
% for the *same* day. 
% Arguably, this could be a between-subjects comparison, looking only at
% speed measurements for days in which divers NEVER entered the water, and 
% comparing to days in which they did. 
% This would make sense if we believed that any time spent with divers 
% alter shark behavior for the remainder of that day. If so, we'd want to
% keep these comparisons close in time (within a day, or few days apart).

dix=  find(Divers == 1);  dix= dix(~isnat( datetimes_sharks(dix) )); % toss NaTs
ndix= find(Divers == 0);  ndix= ndix(~isnat( datetimes_sharks(ndix) )); 

t_Divers= datetimes_sharks(dix);    
t_no_Divers= datetimes_sharks(ndix); 

Speed_d= Speed(dix);    
Speed_nd= Speed(ndix);   

date_Divers= t_Divers;   date_no_Divers= t_no_Divers; 
date_Divers.Format= 'dd-MMM-yyyy'; date_no_Divers.Format= 'dd-MMM-yyyy';

dwd_tmp= [date_Divers.Year date_Divers.Month date_Divers.Day];
dnd_tmp= [date_no_Divers.Year date_no_Divers.Month date_no_Divers.Day];
[dwd_u, inds]= unique(dwd_tmp, 'rows', 'stable'); 

days_w_Divers= date_Divers(inds); % unique days with divers

ndwd= length(days_w_Divers);      % number of days with divers 

tmp_nd_hr= hours(timeofday( t_no_Divers ));
am_feed_win= tmp_nd_hr > 10 & tmp_nd_hr < 11; 
pm_feed_win= tmp_nd_hr > 14.5 & tmp_nd_hr < 15.5;

% specify how you want to run the stats
excludeFeedTimes= 1;
eveningDiversOnly= 0; 

if excludeFeedTimes
    criteria_nd= ~am_feed_win & ~pm_feed_win; 
elseif eveningDiversOnly
    criteria_nd= tmp_nd_hr >= 15.25; 
else % get all times w/o divers from the same day
    criteria_nd= true; 
end

clear spd_dd spd_nd spd_d_and_nd

for dd= 1:ndwd-1
    tmp_dd_dt= t_Divers(inds(dd):inds(dd+1)-1);
    tmp_dd_hr= hours(timeofday(tmp_dd_dt)); 

    tmp_dd_spd= Speed_d(inds(dd):inds(dd+1)); 
    [y, m, d]= ymd( tmp_dd_dt ); 
    same_day= y == y(1) & m == m(1) & d == d(1); 
    spd_dd(dd, 1)= mean( tmp_dd_spd(same_day), 'omitnan' ); 

    sdnd= dnd_tmp(:, 1) == y(1) & dnd_tmp(:, 2) == m(1) & dnd_tmp(:, 3) == d(1);

    thresh_nd= criteria_nd & sdnd; 
    spd_nd(dd, 1)= mean( Speed_nd(thresh_nd), 'omitnan' ); 

end

spd_d_and_nd= [spd_dd spd_nd];           % combine data
ix_nanz= spd_nd == 0 | isnan(spd_nd);    % get rid of missing indices/NaNs
spd_d_and_nd= spd_d_and_nd(~ix_nanz, :); 

summarySTD= std(spd_d_and_nd, 1); 

figure; set(gcf,'color','w');
subplot(1, 2, 1);
if ~exportFig
    rc= withinSubj_raincloud(spd_d_and_nd, rc_col, 'alpha', 1); xlabel(spd_ax_lbl)
else
    rc= withinSubj_raincloud(spd_d_and_nd, rc_col); xlabel(spd_ax_lbl)
end

[h, p, ci, stats]= ttest(spd_d_and_nd(:, 2), spd_d_and_nd(:, 1))

% disp(num2str(p,'%.8f'))

nu = stats.df;
tv = linspace(-15, 15, 300);
tdistpdf = tpdf(tv, nu);
t_obs = stats.tstat;
tvalpdf = tpdf(t_obs, nu);
t_crit = tinv(0.95, nu);

subplot(1, 2, 2); 
plot(tv, tdistpdf); hold on
scatter(t_obs, tvalpdf, "filled")
xline(-t_crit, "--")
legend(["Student's t pdf", "t-Statistic", "Critical Cutoff"])


slopes= spd_d_and_nd(:, 2) - spd_d_and_nd(:, 1); 
pos= slopes > 0;
neg= slopes < 0; 

figure(2); clf; 
ps= plotSpread(slopes); 
[xs, ys]= deal(get(ps{1, 1}, 'XData'), get(ps{1, 1}, 'YData'));

% apply a colormap centered at zero
cs= centerColorsAtZero(ys, m_owt, 0);

figure(2); clf; 
yline(0, '--'); hold on;
scatter(xs, ys, [], ys, 'filled');  colormap(cs); 
xticks(1)
xlim([0 2]);
ylim([-0.4 0.6]);
set(gcf, 'color', 'w');

figure(3); clf; 
circPercent([sum(pos), sum(neg)] ./ length(slopes), 2, 'color', colors); 


%% HELPER FUNCTIONS--------------------------------------------------------

function c= centerColorsAtZero(data, colorMatrix, centerValue)

L=  length(data); 
mx= max(data); 
mi= min(data); 

idx= L * abs(centerValue-mi) / (mx-mi); 

N= size(colorMatrix, 1);

if ~rem(N, 2)         % number of elements is even
    k= (N/2) + 1;
else                  % number of elements is odd
    k= (N+1) / 2; 
end

c1= interp1(1:k, colorMatrix(1:k, :), linspace(1, k, idx), 'linear');
c2= interp1(1:k, colorMatrix(k:N, :), linspace(1, k, L - idx), 'linear');

c=  [c1; c2];

end