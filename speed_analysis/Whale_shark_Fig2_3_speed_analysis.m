%% Whale shark swim speed analysis

% Reproduces panels A-C, & E of Figure 2, panels A-F, of Figure 3, and
%            panels A-B of Figure S5 from supplementary material

% Gallimore C. G.*, Walton C.*, Nugent R., Fradkin M., Poppell L., 
% Schreiber C., Coco C., Grober M., Carlson B., Dove A. D. M., Black M. P. 
% (2024).  A longitudinal behavioral analysis of aquarium whale sharks 
% (Rhincodon typus): insights into individual variation, social hierarchy, 
% and anticipatory cues. Frontiers in Marine Science, 11, 1418002. 
% doi: 10.3389/fmars.2024.1418002)


%% Parameters

computeSpeed= 1;        % make 0 to pull pre-comp'd speeds from csv file
speed_units= 'meters';  % change to feet for ft/sec
pool_morning_evening_for_individuals= 1;   
                        % pool both morning and evening meals for 
                        % breakdown of speed changes by individual 

if strcmpi(speed_units, 'meters')
    spd_ax_lbl= 'Speed (m/sec)'; 
    convF= unitsratio('meters', 'feet');
    speed_bin_edges= linspace(0, 2, 66);
    spd_rho= 0:0.5:2;
    spider_ax= [0.61, 0.13, 0.4, 0.21, 0.04, 0.01; 0.82, 0.25, 0.79, 0.39, 0.19, 0.07];
    spider_ax2= [0.61, 0.13, 0.4, 0.21, 0.04, 0.01, 0.53; 0.82, 0.25, 0.79, 0.39, 0.19, 0.07, 0.98];
    boot_spd_lim= [0.3 1.5];  
    bs_ytix= 0.5:0.25:boot_spd_lim(end);
    bs_ytlbls= num2cell(bs_ytix);
    spd_ci_lim= [-0.15 0.15];
else
    spd_ax_lbl= 'Speed (ft/sec)'; 
    convF= 1;   % keep units same if want in ft/sec
    speed_bin_edges= 0:0.1:6.5;
    spd_rho= 0:2:6; 
    spider_ax= [2, 0.41, 0.4, 0.21, 0.04, 0.01; 2.72, 0.8, 0.79, 0.39, 0.19, 0.07];
    spider_ax2= [2, 0.41, 0.4, 0.21, 0.04, 0.01, 0.53; 2.72, 0.8, 0.79, 0.39, 0.19, 0.07, 0.98];
    boot_spd_lim= [1 5]; 
    bs_ytix= 1:0.5:boot_spd_lim(end);
    bs_ytlbls= num2cell(bs_ytix);
    spd_ci_lim= [-0.45 0.45]; 
end

% make some colormaps
cmap=  cmocean('ice');   

m_owt= flipud( [0.768627450980392	 0.258823529411765	0.101960784313725;
                  0.976470588235294	 0.560784313725490	0.270588235294118;
                  0.996078431372549	 0.749019607843137	0.431372549019608;
                         0.9                0.9                0.9;
                  0.592156862745098	 0.807843137254902	0.800000000000000;
                  0.070588235294117  0.564705882352941	0.556862745098039;
                  0.086274509803922	 0.349019607843137	0.290196078431373] );
    
tealDark= m_owt(2, :);  tealLight= m_owt(3, :); orangeMed= m_owt(6, :); 

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

% missing string values in this dataset input by default as <missing>
% make string for that
miss_str= string( missing );   % find(isnat(datetimes_all))   % get NaT


%% Process datetime data

dates= dataTbl.Date;        missing_dates=   find(ismissing(dates)); 
dates_dt= datetime(dates, 'InputFormat', 'MM/dd/uuuu', 'TimeZone','America/New_York');
TOD= dataTbl.Time_of_Day;   missing_tstamps= find(ismissing(TOD)); 
TOD_dt=   datetime(TOD, 'InputFormat', 'HH:mm');

datetimes_all= dates_dt + timeofday(TOD_dt);  


%% Process shark data and plot some shark descriptives

Shark= upper( dataTbl.Shark );  Shark= Shark(~ismissing(Shark)); % get rid of missing shark

[Shark_g, Shark_Name]= findgroups(Shark);
n_Sharks= length(Shark_Name); 
c_idx= linspace(1, size(cmap, 1), n_Sharks); % define some colors
c_idx2= ceil( linspace(1, size(cmap, 1), n_Sharks+1) ); % define some colors

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

% Process other measures like depth and swim dir
Depth= upper( dataTbl.Depth );  Depth= Depth(~ismissing(Shark));
Dir= upper( dataTbl.Direction_of_Swim );  Dir= Dir(~ismissing(Shark)); 

Feed=  dataTbl.Feed_Status;     

% Initialize categorical variables for depth and swim dir
n_depths= 4;  n_dir= 2; 
d_cats= categories(categorical(Depth));  d_cats= d_cats([3 1 2 4]);
s_dir= categories(categorical(Dir)); 

depth_obs= zeros(n_Sharks, n_depths);  % depth_props= depth_obs; 
swim_dir_obs= zeros(n_Sharks, n_dir);  % swim_dir_props= swim_dir_obs; 

for n= 1:n_Sharks
    % grab speed data
    Shark_logical= strcmpi(Shark_Name(n), Shark);
    Shark_idx{n}=  find(Shark_logical);
    Speed_obs_per_shark{n, 1}= Speed(Shark_idx{n}); 
    Datetimes_per_shark{n, 1}= datetimes_sharks(Shark_idx{n}); 
    dt_hr_each_shark{n, 1}= hours( timeofday( Datetimes_per_shark{n} ) );
    avg_Speed(n)=  mean(Speed_obs_per_shark{n}, 'omitnan');
    std_Speed(n)=  std(Speed(Shark_idx{n}), 'omitnan');

    % grab depth data
    d_sh= reordercats( categorical(Depth(Shark_idx{n})), [3 1 2 4] ); % sort from shallowest to deepest
    sd_sh= categorical(Dir(Shark_idx{n})); 
    nc_d= histcounts(d_sh);  nc_sd= histcounts(sd_sh); 
    depth_obs(n, :)= nc_d;  swim_dir_obs(n, :)= nc_sd(1:n_dir);  
end


%% Plot Figure 2 Panel E

depth_props= depth_obs ./ sum(depth_obs, 2);  % convert to proportions
depth_props(5, :)= nan(1, n_depths);  % pad with nans to make room for lgd

swim_dir_props= swim_dir_obs ./ sum(swim_dir_obs, 2);

cmat= cmap(c_idx(2), :); 

characteristics= [avg_Speed', std_Speed', depth_props(1:4, :), swim_dir_props(:, 2) ];

figure('Name', 'Shark characteristics');    clear sp2
sp2= spider_plot_class(characteristics); 

sp2.AxesPrecision = 2;
sp2.AxesLabels= {'Avg Speed', 'Std Speed', 'Surface', 'Beneath Surface', 'Deep', 'Very Deep', 'CW Pref'}; 
sp2.AxesLimits = spider_ax2;
sp2.Color= cmap(c_idx2(1:4), :);
sp2.LegendLabels = cellstr(Shark_Name)';
sp2.FillOption = 'on';
sp2.FillTransparency =  0.2;


%% Plot Figure S5 Panels A-B

figure('Name', 'Speed histograms all sharks');  

subplot(1, 2, 1);  histogram(Speed, speed_bin_edges, 'FaceAlpha', 1); 
                   xlabel(spd_ax_lbl); ylabel('Count'); xlim([0 2])
                   grid on; set(gca, 'GridAlpha', 1); 
subplot(1, 2, 2);  [hh, sb_mat, x_ax, ~]= stackedHist(Speed, Shark_g, Shark_Name, ...
                                                      'binRange', speed_bin_edges, 'probability'); 
                   xlim([0 2]); colormap(cmap); grid on; set(gca, 'GridAlpha', 1); 


%% Plot Figure 2 Panel A -- Speed CDFs

[spd_ntiles, spd_cdfs]= cellfun(@ecdf, Speed_obs_per_shark, 'UniformOutput', false); 

figure('Name', 'Speed CDFs');   hold on
for n= 1:n_Sharks
    scatter(spd_cdfs{n}, spd_ntiles{n}, 14, cmap(c_idx2(n), :), 'filled');
    if n == n_Sharks
        legend(cellstr(Shark_Name));
        for m= 1:n_Sharks
            [~, loc]= min(abs(avg_Speed(m)- spd_cdfs{m}));
            plot([avg_Speed(m)-std_Speed(m), avg_Speed(m)+std_Speed(m)], ...
                 [spd_ntiles{m}(loc) spd_ntiles{m}(loc)], '-k', 'linewidth', 2, 'HandleVisibility', 'off');
            scatter(avg_Speed(m), spd_ntiles{m}(loc), 60, cmap(c_idx2(m), :), ...
                    'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
        end
    end
end
title('CDF'); grid on
yticks(0:0.25:1); 
ylabel('Quantiles'); xlabel(spd_ax_lbl);
set(gcf, 'Color', 'w'); 

spd_kurt= cellfun(@(x) kurtosis(x, 0), Speed_obs_per_shark, 'UniformOutput', false); 
spd_skew= cellfun(@(x) skewness(x, 0), Speed_obs_per_shark, 'UniformOutput', false); 

% set kurtosis of normal distribution at 0
spd_kurt= cell2mat(spd_kurt) - 3; 


%% Plot Figure 2 Panel B -- Aggregate speed by time of day on a clock

tmin= 0;  tmax= 24;  rotation= 0; 
[tnorm_agg, ctimes]= normalizeTimeCircadianClock(datetimes_sharks, tmin, tmax); 

n_pts= 50; 
tnorm_deg= rad2deg(tnorm_agg);   % 'nbins', [n_pts n_pts]

[~, ~, data, bin_rng]= densityContour(tnorm_deg, Speed, 'nbins', [n_pts n_pts], 'valsOnly', 'probability');  

figure('Name', 'Aggregate speed by time of day'); 
ph= makeClockGrid(30, spd_rho, 'military');

% close the data, ensuring angular vals go all the way 0 to 360, 
% such that data at 360 is the same as data at 0
if logical(mod((bin_rng{1}(1) - bin_rng{1}(end)), 360))
    bin_rng{1}(end+1) = bin_rng{1}(1) + 360;
    data(end+1, :) = data(1, :);    
end

[mt, mspd]= meshgrid(bin_rng{1}, bin_rng{2});
axAng= convStdPolarAngle(mt, 'cw', 'top', rotation);

% Now convert to cartesian coords
[px, py]= pol2cart(axAng*pi/180, mspd); 

% And then use contourf to plot. Optional formatting adjustments follow.
[~, hc]= contourf(px, py, data', 50);
uistack(hc, 'bottom');   % this will put the grid on top of the contours
hc.LevelList(1) = [];    % this removes the contours on lowest level(s) 
hc.LineStyle = 'none';   % enable/disble contour lines
colormap(cmap);  % colormap jet
clim([0 0.008])
title('Speed by time of day: Aggregate')


%% Plot Figure 2 Panel C -- Speed by time of day individually

spd_ctour_all= zeros(n_pts+1, n_pts, n_Sharks); 
spd_bc_all= cell(n_Sharks, 2); 
m_dir_all= permute( spd_ctour_all, [2 1 3] ); 
m_spd_all= m_dir_all; 

px_all= m_dir_all; py_all= m_dir_all; 

figure('Name', 'Individual shark speed by TOD');   
         % for every shark, compute bivariate speed/time density, 
         % "close" angular data
         % place zero (midnight) at the top, and convert to cartesian
for n= 1:n_Sharks
    t_norm_deg_shark{n}= rad2deg( tnorm_agg(Shark_idx{n}) );   % convert to degrees
    t_clock_shark{n}= ctimes(Shark_idx{n});

    % compute density
    [~, ~, spd_ctour, bc_tmp]= densityContour(t_norm_deg_shark{n}, Speed_obs_per_shark{n}, 'nbins', [n_pts n_pts], ...
                                              'valsOnly', 'probability'); 
    % pre-construct clock grid
    subplot(1, n_Sharks, n); 
    ph_ind(n, 1) = makeClockGrid(30, spd_rho, 'military', 'nSubplots', n_Sharks);

    % close the data
    if logical(mod((bc_tmp{1}(1) - bc_tmp{1}(end)), 360))
        bc_tmp{1}(end+1) = bc_tmp{1}(1) + 360;
        spd_ctour(end+1, :) = spd_ctour(1, :);    
    end

    [md_tmp, mspd_tmp] = meshgrid(bc_tmp{1}, bc_tmp{2});

    ang_tmp= convStdPolarAngle(md_tmp, 'cw', 'top', rotation);
    [px_tmp, py_tmp]= pol2cart(ang_tmp*pi/180, mspd_tmp); % convert to cartesian

    [~, hc(n)]= contourf(px_tmp, py_tmp, spd_ctour', 40);
    uistack(hc(n), 'bottom');   % put the grid on top of the contours
    hc(n).LevelList(1) = [];    % remove contour on lowest level(s) 
    hc(n).LineStyle = 'none';   % disable contour lines
    colormap(cmap);  % colormap(cmap3)  % colormap jet % 
    ax_lims{n, 1}= clim; 
    title(strcat('Speed by time of day:', {' '}, Shark_Name(n))); 

    % store data
    spd_ctour_all(:, :, n)= spd_ctour;   spd_bc_all(n, :)= bc_tmp;
    m_dir_all(:, :, n)= md_tmp;          m_spd_all(:, :, n)= mspd_tmp; 
    px_all(:, :, n)= px_tmp;             py_all(:, :, n)= py_tmp; 

end

ax_mins= cellfun(@(x) min(x), ax_lims);  % get color axis lims
ax_maxs= cellfun(@(x) max(x), ax_lims); 

% normalize probability mapping for each subplot
for n= 1:n_Sharks
    subplot(1, 4, n); clim([min(ax_mins) min(ax_maxs)]);
end


%% Process feed times AS A FUNCTION OF SHARKS PRESENT

Feed= Feed(~ismissing(Shark));  

tmp_nf= isnan(Feed); 

% process feed times with NaN before them
nan_onsets= find( diff(tmp_nf) == 1 );  
nan_offsets= find( diff(tmp_nf) == -1 ) +1;  

date_nanON= datetimes_sharks(nan_onsets);
date_nanOFF= datetimes_sharks(nan_offsets); 

n_nan_btw= nan_offsets - nan_onsets; 

time_btw_nan= date_nanOFF - date_nanON; 

% find nan offsets that are a different day to inspect
diffday_nan= date_nanON( hours(time_btw_nan) > 12 );

% process feed times with no NaNs
hungry2fed_ON= find( diff(Feed) == 1 );
hungry2fed_OFF= find( diff(Feed) == 1 ) +1;

time_btw_reg= datetimes_sharks(hungry2fed_OFF)- datetimes_sharks(hungry2fed_ON);

% combine them, sort chronologically, and get feeding duration
all_fON= [nan_onsets; hungry2fed_ON];
all_fOFF= [nan_offsets; hungry2fed_OFF];

[sort_fONs, s_idx]= sort(all_fON); 
sort_fOFFs= all_fOFF(s_idx); 

all_feed_starts= datetimes_sharks(sort_fONs); 

all_feed_durs= datetimes_sharks(sort_fOFFs)- datetimes_sharks(sort_fONs);

% constrain to feed intervals that were less than 1 hour (arbitrary)
sub1hr= hours(all_feed_durs) < 0.75 & minutes(all_feed_durs) > 10;
sub1_ON= sort_fONs(sub1hr); 
sub1_OFF= sort_fOFFs(sub1hr); 
feed_durs_sub1hr= all_feed_durs(sub1hr);
n_feeds= 1:sum(sub1hr); 

% index morning and evening feeds and grab their median interval time
amfs_idx= hours( timeofday(datetimes_sharks(sub1_ON))) < 11;   % morning feed start index
amff_idx= hours( timeofday(datetimes_sharks(sub1_OFF))) < 12;  % morning feed finish index

pmfs_idx= hours( timeofday(datetimes_sharks(sub1_ON))) > 12;   % evening feed start index
pmff_idx= hours( timeofday(datetimes_sharks(sub1_OFF))) > 12;  % evening feed finish index

% for getting pre-meal speed times below
am_feeds_meeting_criteria= sub1_ON(amfs_idx);
pm_feeds_meeting_criteria= sub1_ON(pmfs_idx);
amf_post= sub1_OFF(amff_idx);
pmf_post= sub1_OFF(pmff_idx);

am_feed_starts= datetimes_sharks(am_feeds_meeting_criteria); 
am_feed_fins=   datetimes_sharks(amf_post); 

med_am_feed_int= median(timeofday( [am_feed_starts; am_feed_fins] )); 
med_am_feed_dt=  datetime("today") + med_am_feed_int; % date arbitrary

pm_feed_starts= datetimes_sharks(pm_feeds_meeting_criteria); 
pm_feed_fins=   datetimes_sharks(pmf_post); 

med_pm_feed_int= median(timeofday( [pm_feed_starts; pm_feed_fins] )); 
med_pm_feed_dt=  datetime("today") + med_pm_feed_int; % date arbitrary


%% Analyze speeds leading up to morning meal time
% event-lock the speeds to meal-time and look 90 minutes before

% make time between 9-10:30am the baseline time for morning speed calc
am_tmp_1pad= [0; am_feeds_meeting_criteria]; 

n_am_meals= length(am_tmp_1pad)-1;
n_tpts= 90; 
am_feed_time=  zeros(1, n_am_meals); 
am_90min_b4=   zeros(1, n_am_meals); 
am_resamp_mat= zeros(n_am_meals, n_tpts); 

% get timestamps 90 min before morning feed time for each meal
for am= 1:n_am_meals
    
    am_feed_time(am)= hours( timeofday(datetimes_sharks(am_tmp_1pad(am+1))));
    am_90min_b4(am)= am_feed_time(am) - 1.5; 

    am_resamp_mat(am, :)= linspace(am_90min_b4(am), am_feed_time(am), 90);

    tmp_dt_array=  datetimes_sharks(am_tmp_1pad(am)+1:am_tmp_1pad(am+1)); 
    tmp_spd_array= Speed(am_tmp_1pad(am)+1:am_tmp_1pad(am+1));

    times_up2feed= hours( timeofday(tmp_dt_array) ); 

    [y, m, d]= ymd( tmp_dt_array );
    same_day= y == y(end) & m == m(end) & d == d(end); 

    tmp_hr= times_up2feed > am_90min_b4(am) & same_day;

    am_dt_cell{am, 1}=  tmp_dt_array(tmp_hr)'; 
    am_hr_cell{am, 1}=  hours( timeofday(tmp_dt_array(tmp_hr)) );
    am_spd_cell{am, 1}= tmp_spd_array(tmp_hr);
end

am_resamp_dur= duration(hours(am_resamp_mat), 'format', 'hh:mm');

resamp_sp_am_all= repmat(am_resamp_mat, 1, 1, n_Sharks); 
resamp_t_am_all= resamp_sp_am_all; 

% get speed at resampled time for each shark each meal
for meal= 1:n_am_meals
    ft= am_feed_time(meal);

    for ns= 1:n_Sharks
        sh_ix= Shark_idx{ns}; 
        % get indices within meal bounds
        idx_btw_meal= sh_ix(sh_ix <= am_tmp_1pad(meal+1) & sh_ix >= am_tmp_1pad(meal)+1);
        if isempty(idx_btw_meal)
            resamp_sp_am_all(meal, :, ns)= NaN; 
            resamp_t_am_all(meal, :, ns)= NaN; 
            continue
        else
            sh_sp= Speed(idx_btw_meal);   % shark speed between subsequent meals
            tu2f= hours( timeofday(datetimes_sharks(idx_btw_meal)) ); % times up to feed
            [y, m, d]= ymd( datetimes_sharks(idx_btw_meal));
            same_day= y == y(end) & m == m(end) & d == d(end); 
    
            corr_t= tu2f > am_90min_b4(meal) & same_day; % correct times
    
            t2a= tu2f(corr_t);  % times to analyze
            s2a= sh_sp(corr_t); % speeds to anayze
        end

        for tp= 1:n_tpts-1
            ix_stamp= t2a >= am_resamp_mat(meal, tp) & t2a <= am_resamp_mat(meal, tp+1);

            if any(ix_stamp)
                if sum(ix_stamp) > 1
                    resamp_sp_am_all(meal, tp+1, ns)= mean( s2a(ix_stamp) );
                    resamp_t_am_all(meal, tp+1, ns)= mean( t2a(ix_stamp) );
                else
                    resamp_sp_am_all(meal, tp+1, ns)= s2a(ix_stamp);
                    resamp_t_am_all(meal, tp+1, ns)= t2a(ix_stamp);
                end
            else
                resamp_sp_am_all(meal, tp+1, ns)= NaN; 
                resamp_t_am_all(meal, tp+1, ns)= NaN;
            end
        end
    end
end

resamp_sp_am_all(:, 1, :)= NaN; 
resamp_t_am_all(:, 1, :)= NaN; 

resamp_t_am_idx= mean(resamp_t_am_all, 3, 'omitnan');

rast_rc_am= resamp_t_am_idx; 

tstamp_am_dx= ~isnan(resamp_t_am_idx); 

[r, c]= find(tstamp_am_dx);

for ii= 1:n_am_meals
    real_am_tstamps= c(r == ii)'; 
    resamp_t_am_idx(ii, real_am_tstamps)= real_am_tstamps; 
    rast_rc_am(ii, :)= ii * ones(1, size(rast_rc_am, 2)); 
end

resamp_t_am_col= reshape(resamp_t_am_idx', [], 1); 
rast_col_am= reshape(rast_rc_am', [], 1); 
ts_am_actual= zeros(1, n_am_meals); 
meal_num_am_actual= zeros(1, n_am_meals); 

% grab actual feed time
for f= 1:n_am_meals
    am_actual_ix= find(rast_col_am == f, 1, 'last');
    ts_am_actual(f)= resamp_t_am_col(am_actual_ix); 
    meal_num_am_actual(f)= rast_col_am(am_actual_ix);
end


%% Analyze speeds leading up to afternoon meal time
% event-lock the speeds to meal-time and look 90 minutes before

% make time between 1:30-3pm the baseline time for afternoon speed calc
pm_tmp_1pad= [0; pm_feeds_meeting_criteria]; 

n_pm_meals= length(pm_tmp_1pad)-1;
n_tpts= 90; 
pm_feed_time=  zeros(1, n_pm_meals); 
pm_90min_b4=   zeros(1, n_pm_meals); 
pm_resamp_mat= zeros(n_pm_meals, n_tpts); 

% get timestamps 90 min before afternoon feed time for each meal
for pm= 1:n_pm_meals

    pm_feed_time(pm)= hours( timeofday(datetimes_sharks(pm_tmp_1pad(pm+1))));
    pm_90min_b4(pm)= pm_feed_time(pm) - 1.5; 

    pm_resamp_mat(pm, :)= linspace(pm_90min_b4(pm), pm_feed_time(pm), 90);

    tmp_dt_array=  datetimes_sharks(pm_tmp_1pad(pm)+1:pm_tmp_1pad(pm+1)); 
    tmp_spd_array= Speed(pm_tmp_1pad(pm)+1:pm_tmp_1pad(pm+1));

    times_up2feed= hours( timeofday(tmp_dt_array) ); 

    [y, m, d]= ymd( tmp_dt_array );

    same_day= y == y(end) & m == m(end) & d == d(end); 

    tmp_hr= times_up2feed > pm_90min_b4(pm) & same_day;

    pm_dt_cell{pm, 1}=  tmp_dt_array(tmp_hr)'; 
    pm_hr_cell{pm, 1}=  hours( timeofday(tmp_dt_array(tmp_hr)) );
    pm_spd_cell{pm, 1}= tmp_spd_array(tmp_hr);
end

pm_resamp_dur= duration(hours(pm_resamp_mat), 'format', 'hh:mm');

resamp_sp_pm_all= repmat(pm_resamp_mat, 1, 1, n_Sharks); 
resamp_t_pm_all= resamp_sp_pm_all; 

% get speed at resampled time for each shark each meal
for meal= 1:n_pm_meals
    ft= pm_feed_time(meal);

    for ns= 1:n_Sharks
        sh_ix= Shark_idx{ns}; 
        % get indices within meal bounds
        idx_btw_meal= sh_ix(sh_ix <= pm_tmp_1pad(meal+1) & sh_ix >= pm_tmp_1pad(meal)+1);
        if isempty(idx_btw_meal)
            resamp_sp_pm_all(meal, :, ns)= NaN; 
            resamp_t_pm_all(meal, :, ns)= NaN; 
            continue
        else
            sh_sp= Speed(idx_btw_meal);   % shark speed between subsequent meals
            tu2f= hours( timeofday(datetimes_sharks(idx_btw_meal)) ); % times up to feed
            [y, m, d]= ymd( datetimes_sharks(idx_btw_meal));
            same_day= y == y(end) & m == m(end) & d == d(end); 
    
            corr_t= tu2f > pm_90min_b4(meal) & same_day; % correct times
    
            t2a= tu2f(corr_t);  % times to analyze
            s2a= sh_sp(corr_t); % speeds to anayze
        end

        for tp= 1:n_tpts-1
            ix_stamp= t2a >= pm_resamp_mat(meal, tp) & t2a <= pm_resamp_mat(meal, tp+1);

            if any(ix_stamp)
                if sum(ix_stamp) > 1
                    resamp_sp_pm_all(meal, tp+1, ns)= mean( s2a(ix_stamp) );
                    resamp_t_pm_all(meal, tp+1, ns)= mean( t2a(ix_stamp) );
                else
                    resamp_sp_pm_all(meal, tp+1, ns)= s2a(ix_stamp);
                    resamp_t_pm_all(meal, tp+1, ns)= t2a(ix_stamp);
                end
            else
                resamp_sp_pm_all(meal, tp+1, ns)= NaN; 
                resamp_t_pm_all(meal, tp+1, ns)= NaN;
            end
        end
    end
end

resamp_sp_pm_all(:, 1, :)= NaN; 
resamp_t_pm_all(:, 1, :)= NaN; 

resamp_t_pm_idx= mean(resamp_t_pm_all, 3, 'omitnan');

rast_rc_pm= resamp_t_pm_idx; 

tstamp_pm_dx= ~isnan(resamp_t_pm_idx); 

[r, c]= find(tstamp_pm_dx);

for ii= 1:n_pm_meals
    real_pm_tstamps= c(r == ii)'; 
    resamp_t_pm_idx(ii, real_pm_tstamps)= real_pm_tstamps; 
    rast_rc_pm(ii, :)= ii * ones(1, size(rast_rc_pm, 2)); 
end

resamp_t_pm_col= reshape(resamp_t_pm_idx', [], 1); 
rast_col_pm= reshape(rast_rc_pm', [], 1); 
ts_pm_actual= zeros(1, n_pm_meals); 
meal_num_pm_actual= zeros(1, n_pm_meals); 

% grab actual feed time
for f= 1:n_pm_meals
    pm_actual_ix= find(rast_col_pm == f, 1, 'last');
    ts_pm_actual(f)= resamp_t_pm_col(pm_actual_ix); 
    meal_num_pm_actual(f)= rast_col_pm(pm_actual_ix);
end


%% Plot Figure 3 Panel F -- Speed change by individual animal for am & pm

b_ix= 11:20;  % base window
f_ix= 81:90;  % feed window

% compute baseline and feed for morning & evening matrices
ind_sh_base_am= squeeze(mean(resamp_sp_am_all(:, b_ix, :), 2, 'omitnan'));
ind_sh_feed_am= squeeze(mean(resamp_sp_am_all(:, f_ix, :), 2, 'omitnan'));
ind_sh_base_pm= squeeze(mean(resamp_sp_pm_all(:, b_ix, :), 2, 'omitnan'));
ind_sh_feed_pm= squeeze(mean(resamp_sp_pm_all(:, f_ix, :), 2, 'omitnan'));

% calculate mean for each shark in both time windows
ind_mean_am= [mean(ind_sh_base_am, 1, 'omitnan'); mean(ind_sh_feed_am, 1, 'omitnan')];
ind_mean_pm= [mean(ind_sh_base_pm, 1, 'omitnan'); mean(ind_sh_feed_pm, 1, 'omitnan')];

% create baseline & feed x coordinates
xb_am= repmat(1:2:(n_Sharks*2)-1, n_am_meals, 1);
xf_am= repmat(2:2:n_Sharks*2, n_am_meals, 1);
xb_pm= repmat(1:2:(n_Sharks*2)-1, n_pm_meals, 1);
xf_pm= repmat(2:2:n_Sharks*2, n_pm_meals, 1);

% concatenate coordinates
xc_am= [xb_am(:)'; xf_am(:)'];
xc_pm= [xb_pm(:)'; xf_pm(:)'];
yc_am= [ind_sh_base_am(:)'; ind_sh_feed_am(:)']; 
yc_pm= [ind_sh_base_pm(:)'; ind_sh_feed_pm(:)']; 

% compute change in speed
slopes_am= ind_sh_feed_am - ind_sh_base_am; 
slopes_pm= ind_sh_feed_pm - ind_sh_base_pm; 

% store total meal counts for each shark
n_meals_scored_am= sum(~isnan(slopes_am), 1); 
n_meals_scored_pm= sum(~isnan(slopes_pm), 1); 

% mean line coords
xm_am= reshape(1:n_Sharks*2, 2, n_Sharks); 
xm_pm= xm_am; 

% plot everything
if pool_morning_evening_for_individuals == 1  % if pooling data

    xc= [xc_am, xc_pm];  yc= [yc_am, yc_pm];  xm= xm_am; 

    ind_mean= [mean([ind_sh_base_am; ind_sh_base_pm], 1, 'omitnan'); ...
               mean([ind_sh_feed_am; ind_sh_feed_pm], 1, 'omitnan')];

    slopes= [slopes_am; slopes_pm]; 

    n_meals_scored= n_meals_scored_am + n_meals_scored_pm; 

    figure('Name', 'Pooled meal speeds by individual'); 
    set(gcf, 'color', 'w')
    s1= subplot(4, 1, 1:2); hold on; 
    plot(xc, yc, '-', 'color', [0.5 0.5 0.5]);
    plot(xm, ind_mean, '-k', 'LineWidth', 2);
    ylabel(spd_ax_lbl);
    xlim([0 n_Sharks*2+1])
    
    s2= subplot(4, 1, 3:4); hold on; 
    ps= plotSpread(slopes, 'xValues', 1.5:2:7.5);
    sc= [get(ps{1, 1}, 'XData'); get(ps{1, 1}, 'YData')];
    
    xs_tmp= sc(1:4);   xs= [xs_tmp{:}];
    ys_tmp= sc(5:8);   ys= [ys_tmp{:}];
    
    % apply a colormap centered at zero
    cs= centerColorsAtZero(ys, m_owt, 0); 
    
    cla(s2); 
    s2= subplot(4, 1, 3:4); hold on; 
    yline(0, '--'); 
    scatter(xs, ys, [], ys, 'filled');  colormap(cs); 
    addSummaryStatLines(slopes, 4, 'median', 'xVals', 1.5:2:7.5);
    ylabel(strcat('\Delta', spd_ax_lbl))
    xlim([0 n_Sharks*2+1])
%     colorbar
    
    props_spd_inc= [sum(slopes > 0, 1) ./ n_meals_scored; 
                    1-(sum(slopes > 0, 1) ./ n_meals_scored)]; 
    
else
    f1= figure('Name', 'AM meal speeds by individual');   set(gcf, 'color', 'w')
    f2= figure('Name', 'PM meal speeds by individual');   set(gcf, 'color', 'w')

    xc= {xc_am, xc_pm};  yc= {yc_am, yc_pm};  xm= xm_am; 

    ind_mean= {ind_mean_am; ind_mean_pm};

    slopes= {slopes_am; slopes_pm}; 

    n_meals_scored= [n_meals_scored_am; n_meals_scored_pm]; 

    for m= 1:2
        if m == 1
            figure(f1); 
        else
            figure(f2);
        end

        s1(m)= subplot(4, 1, 1:2); hold on; 
        plot(xc{m}, yc{m}, '-', 'color', [0.5 0.5 0.5]);
        plot(xm, ind_mean{m}, '-k', 'LineWidth', 2);
        ylabel(spd_ax_lbl);
        xlim([0 n_Sharks*2+1])
        
        s2(m)= subplot(4, 1, 3:4); hold on; 
        ps= plotSpread(slopes{m}, 'xValues', 1.5:2:7.5);
        sc= [get(ps{1, 1}, 'XData'); get(ps{1, 1}, 'YData')];
        
        xs_tmp= sc(1:4);   xs= [xs_tmp{:}];
        ys_tmp= sc(5:8);   ys= [ys_tmp{:}];
        
        % apply a colormap centered at zero
        cs= centerColorsAtZero(ys, m_owt, 0); 
        
        cla(s2(m)); 
        s2(m)= subplot(4, 1, 3:4); hold on; 
        yline(0, '--'); 
        scatter(xs, ys, [], ys, 'filled');  colormap(cs); 
        addSummaryStatLines(slopes{m}, 4, 'median', 'xVals', 1.5:2:7.5);
        ylabel(strcat('\Delta', spd_ax_lbl))
        xlim([0 n_Sharks*2+1])
        
        props_spd_inc= [sum(slopes{m} > 0, 1) ./ n_meals_scored(m, :); 
                        1-(sum(slopes{m} > 0, 1) ./ n_meals_scored(m, :))]; 
         

        compare_am_pm(:, m)=  sum([sum(slopes{m} > 0, 1); n_meals_scored(m, :)], 2);
    end

    [h, p, stats]= fishertest(compare_am_pm);

end


%% Number of feeds am vs pm stats

obs_per_meal_am= sum( sum( ~isnan(resamp_t_am_all), 3, 'omitnan'), 2);
obs_per_meal_pm= sum( sum( ~isnan(resamp_t_pm_all), 3, 'omitnan'), 2);

c_obs= {obs_per_meal_am, obs_per_meal_pm};
c_grp= [1, 2];
def_c= [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
mrkrs= ["square", "^"];

figure('Name', 'N obs per feed, morning vs evening');
subplot(1, 2, 1)
histogram(obs_per_meal_am, 'BinMethod', 'sturges'); hold on;
histogram(obs_per_meal_pm, 'BinMethod', 'sturges');
subplot(1, 2, 2); hold on
po= plotSpread(c_obs, 'distributionMarkers', mrkrs, 'distributionColors', def_c, 'xValues', c_grp);
addSummaryStatLines(c_obs, 2, 'mean');

[h, p, ci, stats]= ttest2(obs_per_meal_am, obs_per_meal_pm);


%% Get final means and linear changepoints

ws_am_mean= mean( resamp_sp_am_all(:, 2:end, :), 3, 'omitnan' ); 
ws_pm_mean= mean( resamp_sp_pm_all(:, 2:end, :), 3, 'omitnan' ); 

% raw traces
% figure;
% subplot(2, 1, 1);
% plot(ws_am_mean');
% subplot(2, 1, 2);
% plot(ws_pm_mean');

an2= sum(isnan(ws_am_mean), 1);  an3= sum(isnan(ws_pm_mean), 1); 

tp_am_mean= mean(ws_am_mean, 1, 'omitnan'); 
tp_pm_mean= mean(ws_pm_mean, 1, 'omitnan'); 
tp_am_SEM= std(ws_am_mean, 1, 'omitnan') ./ sqrt(n_am_meals-an2); 
tp_pm_SEM= std(ws_pm_mean, 1, 'omitnan') ./ sqrt(n_pm_meals-an3); 

smth_am_mean= smooth(tp_am_mean, 10); 
smth_am_SEM=  smooth(tp_am_SEM, 10); 
smth_pm_mean= smooth(tp_pm_mean, 10); 
smth_pm_SEM=  smooth(tp_pm_SEM, 10); 

% compute linear changepoints
[ipt_am_lin, res_am_lin] = findchangepts(tp_am_mean, 'MaxNumChanges', 1, 'Statistic', 'linear');
[ipt_am_smthlin, res_am_smthlin] = findchangepts(smth_am_mean, 'MaxNumChanges', 1, 'Statistic', 'linear');

[ipt_pm_lin, res_pm_lin] = findchangepts(tp_pm_mean, 'MaxNumChanges', 1, 'Statistic', 'linear');
[ipt_pm_smthlin, res_pm_smthlin] = findchangepts(smth_pm_mean, 'MaxNumChanges', 1, 'Statistic', 'linear');

shadeTOP= max([n_am_meals n_pm_meals]) + 10;
xp = [90, 90, 100, 100];
yp = [-10, shadeTOP, shadeTOP, -10];     


%% Plot Figure 3 Panels A-D

figure;  set(gcf,'color','w');
subplot(3, 2, 1); histogram(feed_durs_sub1hr, 35, 'FaceAlpha', 1); 
                  grid on;  set(gca, 'GridAlpha', 1); 
                  ylabel('Count'); xlabel('Feed duration (min)');
                  title('quantified by elapsed time before & after')
subplot(3, 2, 2); makeClockScatter(datetimes_sharks(sub1_ON), n_feeds, 0, 24, 12, ...
                                   'DotSize', 10, 'Color', 'g', 'EdgeColor', 'none', ...
                                   'inclHand', med_am_feed_dt); 
                  set(gca, 'GridAlpha', 1);  hold on; 
                  makeClockScatter(datetimes_sharks(sub1_OFF), n_feeds, 0, 24, 12, ...
                                   'DotSize', 10, 'Color', 'r', 'EdgeColor', 'none', ...
                                   'inclHand', med_pm_feed_dt); 
                  set(gca, 'GridAlpha', 1);
                  title('feed \color{green}onsets \color{black}and \color{red}offsets')
subplot(3, 2, 3); plot_raster_event( resamp_t_am_col, rast_col_am );  hold on; 
                  plot_raster_event( ts_am_actual, meal_num_am_actual, [], '^r' ); 
                  fill(xp, yp, [0 1 0]);
                  xlim([0 95]);  ylabel('meals'); ylim([0 125])
                  title('timestamps leading up to feed time (am)')
subplot(3, 2, 4); plot_raster_event( resamp_t_pm_col, rast_col_pm );  hold on; 
                  plot_raster_event( ts_pm_actual, meal_num_pm_actual, [], '^r' ); 
                  fill(xp, yp, [0 1 0]);
                  xlim([0 95]);  ylabel('meals'); ylim([0 100])
                  title('timestamps leading up to feed time (pm)')
subplot(3, 2, 5); shadedErrorBar([], tp_am_mean, tp_am_SEM, 'transparent', 0);  hold on;
                  shadedErrorBar([], smth_am_mean, smth_am_SEM, 'lineprops','-b', 'transparent', 0); 
                  plot(ipt_am_smthlin, smth_am_mean(ipt_am_smthlin), '^r', 'MarkerSize', 8, 'LineWidth', 3, 'MarkerEdgeColor', 'r');  
                  ylabel(spd_ax_lbl); xlabel('timestamp (min)'); xlim([0 95]);
                  title('event-locked speed before morning meal')
                  legend(["time-pt means", "smoothed (10 min MA)", "change point (linear)"], 'Location', 'northwest'); 
subplot(3, 2, 6); shadedErrorBar([], tp_pm_mean, tp_pm_SEM, 'transparent', 0);  hold on;
                  shadedErrorBar([], smth_pm_mean, smth_pm_SEM, 'lineprops','-b', 'transparent', 0); 
                  plot(ipt_pm_smthlin, smth_pm_mean(ipt_pm_smthlin), '^r', 'MarkerSize', 8, 'LineWidth', 3, 'MarkerEdgeColor', 'r'); 
                  ylabel(spd_ax_lbl); xlabel('timestamp (min)'); xlim([0 95]);
                  title('event-locked speed before evening meal')
                  legend(["time-pt means", "smoothed (10 min MA)", "change point (linear)"], 'Location', 'northwest'); 

% hgexport(gcf, 'meal_speed.eps')


%% Plot Figure 3 Panel E -- After bootstrapping big 3D matrix

colors= [orangeMed; tealDark]; 
nreps= 3000; 
ran_seed= 1; 

baseline_am= resamp_sp_am_all(:, 11:20, :);
beforefd_am= resamp_sp_am_all(:, 81:end, :);

baseline_pm= resamp_sp_pm_all(:, 11:20, :);
beforefd_pm= resamp_sp_pm_all(:, 81:end, :);

npts_base_am= min( min(sum(~isnan(baseline_am), 1), [], 3) ); 
npts_feed_am= min( min(sum(~isnan(beforefd_am), 1), [], 3) ); 

npts_base_pm= min( min(sum(~isnan(baseline_pm), 1), [], 3) ); 
npts_feed_pm= min( min(sum(~isnan(beforefd_pm), 1), [], 3) ); 

npts_boot_am= min([npts_base_am npts_feed_am]);
npts_boot_pm= min([npts_base_pm npts_feed_pm]);

% bootstrap resampled matrices for each shark over specified windows
[bootbase_mat_am, bootfeed_mat_am] = bootstrap_same_indices(baseline_am, beforefd_am, nreps, n_Sharks);
[bootbase_mat_pm, bootfeed_mat_pm] = bootstrap_same_indices(baseline_pm, beforefd_pm, nreps, n_Sharks);

% average over matrices that were sampled with the same indices
bootbaseAV_am= mean( bootbase_mat_am, 3, 'omitnan' ); 
bootfeedAV_am= mean( bootfeed_mat_am, 3, 'omitnan' ); 
bootbaseAV_pm= mean( bootbase_mat_pm, 3, 'omitnan' ); 
bootfeedAV_pm= mean( bootfeed_mat_pm, 3, 'omitnan' ); 

% now randomly sample sharks to include in aggregate average
rng(ran_seed, 'twister'); 
r_subj= randi(n_Sharks, nreps, n_Sharks);

boot_Bshuff_am=  zeros(n_Sharks, nreps);  boot_Fshuff_am=  zeros(n_Sharks, nreps); 
boot_Bshuff_pm=  zeros(n_Sharks, nreps);  boot_Fshuff_pm=  zeros(n_Sharks, nreps); 

for rb= 1:nreps
    boot_Bshuff_am(:, rb)= bootbaseAV_am(r_subj(rb, :), rb); 
    boot_Fshuff_am(:, rb)= bootfeedAV_am(r_subj(rb, :), rb); 
    boot_Bshuff_pm(:, rb)= bootbaseAV_pm(r_subj(rb, :), rb); 
    boot_Fshuff_pm(:, rb)= bootfeedAV_pm(r_subj(rb, :), rb);
end

boot_Bshuff_am= boot_Bshuff_am';   boot_Fshuff_am= boot_Fshuff_am';
boot_Bshuff_pm= boot_Bshuff_pm';   boot_Fshuff_pm= boot_Fshuff_pm';

% average over sharks
boot_muB_am= mean(boot_Bshuff_am, 2);  boot_muF_am= mean(boot_Fshuff_am, 2); 
boot_muB_pm= mean(boot_Bshuff_pm, 2);  boot_muF_pm= mean(boot_Fshuff_pm, 2); 

% compute probabilities
[prob1, p_mat1, p_ctrs1]= get_direct_prob(boot_muF_am, boot_muB_am); 
[prob2, p_mat2, p_ctrs2]= get_direct_prob(boot_muF_pm, boot_muB_pm); 

boot_comps_all= [boot_muB_am boot_muF_am boot_muB_pm boot_muF_pm]; 
colors2= flipud( [tealDark; [.5 .5 .5]; tealDark; [.5 .5 .5]] ); 

figure('Name', 'Bootstrapped distributions'); set(gcf,'color','w');
for n= 1:4
    hv{n}=  Violin(boot_comps_all(:, n), n, 'violinColor', colors2(n, :), 'boxcolor', [0 0 0], 'boxwidth', 0.025); hold on
    if n == 2 || n == 4
        xC= [get(hv{n-1}.MedianPlot, 'XData'); get(hv{n}.MedianPlot, 'XData')];
        yC= [get(hv{n-1}.MedianPlot, 'YData'); get(hv{n}.MedianPlot, 'YData')];
        hv{2, n}= plot(xC, yC, '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5);  
    end
end
title(['boostrapped distributions baseline vs feed, N=', num2str(nreps), ' ', 'iters'])
xlim([0.5 4.5]); ylim(boot_spd_lim)
yticks(bs_ytix); yticklabels(bs_ytlbls)
xticks(1:4); xticklabels({'am base', 'am feed', 'pm base', 'pm feed'})
ylabel(spd_ax_lbl)
grid on

boot_means= mean(boot_comps_all, 1);
boot_SEM= std(boot_comps_all, 1);


%% HELPER FUNCTIONS--------------------------------------------------------

function c= centerColorsAtZero(data, colorMatrix, centerValue)

L=  length(data); 
mx= max(data); 
mi= min(data); 

idx= L * abs(centerValue-mi) / (mx-mi); 

N= size(colorMatrix, 1);

if ~rem(N, 2)         % number of elements is even
    k=  N/2;
    k2= k+1;
else                  % number of elements is odd
    k= (N+1) / 2; 
    k2= k; 
end

c1= interp1(1:k, colorMatrix(1:k, :), linspace(1, k, idx), 'linear');
c2= interp1(1:k, colorMatrix(k2:N, :), linspace(1, k, L - idx), 'linear');

c=  [c1; c2];

end