%% Whale shark depth analysis


% make some colormaps
cmap= cmocean('ice');   
cmap2= cmocean('matter');  cm= flipud(cmap2); 


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

Shark= upper( dataTbl.Shark );  Shark= Shark(~ismissing(Shark)); % get rid of missing shark

[Shark_g, Shark_Name]= findgroups(Shark);
n_Sharks= length(Shark_Name); 

datetimes_sharks= datetimes_all(~ismissing(Shark));  % get all datetimes where shark was noted


% Process other measures like depth and swim dir
Depth= upper( dataTbl.Depth );  Depth= Depth(~ismissing(Shark));

for n= 1:n_Sharks
    % grab some data
    Shark_logical= strcmpi(Shark_Name(n), Shark);
    Shark_idx{n}=  find(Shark_logical);
    Datetimes_per_shark{n, 1}= datetimes_sharks(Shark_idx{n}); 
    dt_hr_each_shark{n, 1}= hours( timeofday( Datetimes_per_shark{n} ) );
end


% Plot depth for each shark and number of observations
n_depths= 4;  n_dir= 2; 
d_cats= categories(categorical(Depth));  d_cats= d_cats([3 1 2 4]);

depth_obs= zeros(n_Sharks, n_depths);  % depth_props= depth_obs; 

for n= 1:n_Sharks
    idx= Shark_idx{n};  nObs= length(idx); 
    d_sh= reordercats( categorical(Depth(Shark_idx{n})), [3 1 2 4] ); % sort from shallowest to deepest
    nc_d= histcounts(d_sh);  
    depth_obs(n, :)= nc_d;  
%     depth_props(n, :)= nc_d ./ sum(nc_d);  swim_dir_props(n, :)= nc_sd; 
end

depth_props= depth_obs ./ sum(depth_obs, 2);  % convert to proportions
depth_props(5, :)= nan(1, n_depths);  % pad with nans to make room for lgd

[yt, xt]= ndgrid(1:n_depths, 1:n_depths); 
diams= flipud( depth_props(1:4, :)' * 100 ); 

c_scale= round( max_normalize(diams(:)')*255 ) + 1;

figure('Name', 'Depth props each shark 2'); 
bubblechart(xt(:), yt(:), diams(:), cm(c_scale, :), 'MarkerFaceAlpha', 1); 
bubblesize([10 60])
xlim([0 5]); ylim([0.5 4.5]);
xticks(1:4); xticklabels(Shark_Name)
yticks(1:4); yticklabels(flipud(d_cats))
bubblelegend('N Obs at Depth','Location','eastoutside')
set(gca, 'XAxisLocation', 'top'); set(gcf,'color','w');


%% Omnibus chi-square depth

% depth_obs= depth_obs'; 

rt= sum(depth_obs, 2);     % sum rows
ct= sum(depth_obs, 1);     % sum cols 
d_totals= sum(rt);         % sum total

d_exp= rt * ct / d_totals; 

d_stats= X2ind(depth_obs, d_exp, 0.05); 

D_stat_mat= d_stats.Computed;

omni_p_vals= chi2cdf( D_stat_mat, 9, 'upper'); 
pd_chi_agg=  chi2cdf( sum(d_stats.X2), 9, 'upper'); 

% penalize p_thresh by number of components in omnibus chi-square
omni_p_thresh= 0.05/numel(depth_obs);

% random numbers from the below chi distribution will only be greater than
% 'omni_X2_thresh' 5% of the time
omni_X2_thresh= chi2inv(1-omni_p_thresh, 9);

% zero non-significant X2 components (nan doesnt work for this case)
d_stat_dir= D_stat_mat;
d_stat_dir(d_stat_dir < omni_X2_thresh)= 0;

% make values where observed < expected negative to show directionality
d_stat_dir(D_stat_mat > omni_X2_thresh & depth_obs < d_exp)= -d_stat_dir(D_stat_mat > omni_X2_thresh & depth_obs < d_exp); 

% uncomment to plot p-vals
% figure('Name', 'Omnibus chi-square'); 
% pd_table2= stat_heatmap_directional(d_stat_dir, omni_X2_thresh, 'rowlabels', Shark_Name, 'collabels', d_cats, ...
%                                     'XaxisLocation', 'top', 'pvals', {omni_p_vals, omni_p_thresh});

% plot X2 statistics
figure('Name', 'Omnibus chi-square'); 
pd_table2= statHeatmapDirectional(d_stat_dir, omni_X2_thresh, 'rowlabels', Shark_Name, ...
                                  'collabels', d_cats, 'XaxisLocation', 'top');
colormap(redblue)
pd_table2.colorbar.Label.String= 'X^{2} magnitude'; 
pd_table2.colorbar.Label.FontWeight= 'bold'; 
clim([-600 600])
xlabel('Depth');  ylabel('Shark');
set(gcf,'color','w');


%% Pairwise analysis

pw_comp= nchoosek(1:n_Sharks, 2);
d_comp= nchoosek(1:n_depths, 2);

n_poss_shark_comps= size(pw_comp, 1);
n_poss_depth_comps= size(d_comp, 1);

pwp_all= n_poss_depth_comps * n_poss_shark_comps; 

depth_obs_d= depth_obs;

m_DDNf= zeros(n_Sharks*n_depths); m_DDNf2= m_DDNf; m_DDNf3= m_DDNf; 
D_pw_stat_all= m_DDNf;  D_pw_df_all= m_DDNf;  m_DDN= m_DDNf;


% pairwise depth comparisons 1 shark vs all others
for d_s= 1:n_depths

    d2c= ismember(1:n_depths, d_s);  % depths to compare
    
    for comp_s= 1:size(pw_comp, 1)
        
        s2c= depth_obs_d(pw_comp(comp_s, :), :); % sharks to compare
        
        d2c_dep= s2c(:, d2c);
        comp_dep= sum(s2c(:, ~d2c), 2);

        md= [d2c_dep comp_dep]; 

        dc_mat_all= md;
        pwr_all= sum(dc_mat_all, 2);     % sum rows
        pwc_all= sum(dc_mat_all, 1);     % sum cols 
        pw_totals_all= sum(pwr_all);     % sum total
        
        dpw_exp_all= ceil( pwr_all * pwc_all / pw_totals_all ); 
        
        dpw_stats_all= X2ind(dc_mat_all, dpw_exp_all, 0.05); 
        
        D_pw_stat_all(pw_comp(comp_s, 1) + ((d_s-1)*4), pw_comp(comp_s, 2) + ((d_s-1)*4))= sum( dpw_stats_all.Computed(:) );
        D_pw_df_all(pw_comp(comp_s, 1) + ((d_s-1)*4), pw_comp(comp_s, 2) + ((d_s-1)*4))= dpw_stats_all.df;

        pv= chi2cdf( D_pw_stat_all(pw_comp(comp_s, 1) + ((d_s-1)*4), pw_comp(comp_s, 2) + ((d_s-1)*4)), dpw_stats_all.df, 'upper');
        
        % store chi2cdf p-values
        m_DDN(pw_comp(comp_s, 1) + ((d_s-1)*4), pw_comp(comp_s, 2) + ((d_s-1)*4))= pv; 

        
        [dec, post_p, fstats]= fishertest(dc_mat_all, 'Tail', 'right', 'Alpha', 0.05/(n_Sharks*n_depths) ); 
        m_DDNf(pw_comp(comp_s, 1) + ((d_s-1)*4), pw_comp(comp_s, 2) + ((d_s-1)*4))= pv; 

        [dec, post_p, fstats]= fishertest(dc_mat_all, 'Tail', 'left', 'Alpha', 0.05/(n_Sharks*n_depths) ); 
        m_DDNf2(pw_comp(comp_s, 1) + ((d_s-1)*4), pw_comp(comp_s, 2) + ((d_s-1)*4))= post_p; 
% 
        [dec, post_p, fstats]= fishertest(dc_mat_all, 'Alpha', 0.05/(n_Sharks*n_depths) ); 
        m_DDNf3(pw_comp(comp_s, 1) + ((d_s-1)*4), pw_comp(comp_s, 2) + ((d_s-1)*4))= post_p; 
    
    end
end

p_thresh= (0.05/(2*2))/(size(pw_comp, 1)*n_depths);
X2_thresh= chi2inv(1-p_thresh, 1);


d_stat_dir2= D_pw_stat_all;
d_stat_dir2(D_pw_stat_all > X2_thresh & m_DDNf2 < p_thresh)= -d_stat_dir2(D_pw_stat_all > X2_thresh & m_DDNf2 < p_thresh); 
d_stat_dir2(m_DDNf2 < m_DDNf & d_stat_dir2 > 0)= -d_stat_dir2(m_DDNf2 < m_DDNf & d_stat_dir2 > 0); 

col_lbls= repmat(Shark_Name', 1, n_depths);
row_lbls= repmat(Shark_Name, n_depths, 1);

m_DDNf(m_DDNf == 0)= 1; 


% plot pairwise comparisons -- dotplot
ut= triu(d_stat_dir2, 1);
ut(~triu(d_stat_dir2, 1))= nan; 
ut_p= chi2cdf( abs(ut), dpw_stats_all.df, 'upper'); 


figure('Name', 'Pairwise comparisons'); 
h=  xcorrDotPlot(d_stat_dir2, 'upper', 1, 'rowlabels', row_lbls, 'collabels', col_lbls, ...
                 'majorgrid', 0, 'minorgrid', 1); 
colormap(redblue)
clim([-600 600])


%% Depth occupation diversity

n_obs= sum(depth_obs, 2);
ci= ceil( linspace(1, size(cmap, 1), n_Sharks+1) ); % define some colors

% shannon's diversity index
% higher uncertainty -> higher diversity
% increasing uniformity of the distribution -> increases diversity
d_p= depth_props(1:4, :); 
SDI= shannonEntropy(d_p);

% SDI is equivalent to
% ((n_obs .* log2(n_obs)) - sum(depth_obs .* log2(depth_obs), 2)) ./ n_obs
% when using frequency counts

H_max= log2(n_depths);

% homogeneity / shannon equitability index
HG= SDI ./ H_max; 

HG2= (SDI - min(SDI)) ./ (H_max - min(SDI));

% simpson's similarity index (D)
% probability of randomly picking 2 behaviors belonging to the same
% category
D= sum(depth_obs .* (depth_obs-1), 2) ./ (n_obs .* (n_obs-1));

% plot diversity/similarity indices
figure; 
subplot(2, 1, 1)
yl= yline(H_max, '--k', 'max diversity'); hold on;
scatter(1:n_Sharks, SDI, 200, 1:n_Sharks, 'filled'); 

% set box property to off and remove background color
ax= gca; 
set(ax, 'box', 'off', 'color', 'none', 'TickDir', 'out')
ax.XAxis.Visible= 'off';
yl.LabelHorizontalAlignment = 'right';
yl.LabelVerticalAlignment = 'bottom';

xlim([0 n_Sharks+1])
ylim([1 2])
ylabel('shannon diversity')
grid on; 

subplot(2, 1, 2)
scatter(1:n_Sharks, D, 200, 1:n_Sharks, 'filled', 'Marker', 'diamond');
colormap(cmap(ci(1:4), :))

ax= gca; 
set(ax, 'box', 'off', 'color', 'none', 'TickDir', 'out')
ax.XAxis.Visible= 'off';

xlim([0 n_Sharks+1])
ylim([0 1])
ylabel('simpson similarity')
grid on; 
set(gcf, 'color', 'w'); 
