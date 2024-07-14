%% Whale shark lead-follow interactions


% relations key:
% A = Ahead
% N = Next to
% B = Behind
% X = Crossed

% depth relations
% L = Lower
% E = Equivalent
% H = Higher
% D = Very deep


% make some colormaps
cmap=  cmocean('tempo');  
cmap2= cmocean('ice');  

% parameters
useComputed= 1; 
writeProcTable= 0; 
saveDir= pwd;


%% Load data

fileName= 'GAA_whale_shark_master_datasheet.csv'; 

opts= detectImportOptions( fileName );

opts= setvartype(opts, {'Day', 'Observer', 'Location', ...
                        'Shark', 'Depth', 'Direction_of_Swim', ...
                        'Date', 'Time_of_Day', 'Away_from_Wall', ...
                        'Raw_w_another_shark', 'Raw_Position_other_shark', ...
                        'Raw_Depth_other_shark', 'Other_Notes'}, 'string');      

dataTbl = readtable(fileName, opts);   

Shark= upper( dataTbl.Shark );  Shark= Shark(~ismissing(Shark)); % get rid of missing shark


%% Make a processed dataTbl to manually standardize data key

% grab all necessary cells of the main datasheet to process L-F interxns
dataLF= dataTbl(:, ["Shark", "Direction_of_Swim", "Depth", ...
                    "Raw_w_another_shark", "Raw_Position_other_shark", "Raw_Depth_other_shark", ...
                    "Second_other_Shark", "Second_Shark_Position", "Second_Shark_Depth", ...
                    "Third_Other_Shark", "Third_Shark_Position", "Third_Shark_Depth"]);

other_shark= dataLF.Raw_w_another_shark; 

a= dataTbl.Raw_w_another_shark;

empt_shark= ismissing(a);   

% the below indices correspond to the rows in the main datasheet containing
% lead-follow interactions for manual inspection by authors MB and CW
LF_dataIndices= find(~ismissing(other_shark) & ~empt_shark); 

procLF= dataLF(LF_dataIndices, :); 
if writeProcTable == 1
    fname= strcat(char(datetime("today", 'Format', 'yyMMdd')), '_LF_to_change.xls');  
    fsave= fullfile(saveDir, fname);
    writetable(procLF, fsave, 'Sheet', 1);
end


%% Load cleaned qualitative lead-follow data

fn= '240112_LeadFollow_MBnCW.csv'; 

opts= detectImportOptions( fn );

opts= setvartype(opts, {'Shark', 'Depth', 'Direction_of_Swim', ...
                        'Raw_w_another_shark', 'Raw_Position_other_shark', 'Raw_Depth_other_shark', ...
                        'Second_other_Shark', 'Second_Shark_Position', 'Second_Shark_Depth', ...
                        'Third_Other_Shark', 'Third_Shark_Position', 'Third_Shark_Depth'}, 'string');      

procTbl= readtable(fn, opts);   

Shark= upper( procTbl.Shark );
[Shark_g, Shark_Name]= findgroups(Shark);

other_sharks= [procTbl.Raw_w_another_shark, procTbl.Second_other_Shark, procTbl.Third_Other_Shark]; 
positions= [procTbl.Raw_Position_other_shark, procTbl.Second_Shark_Position, procTbl.Third_Shark_Position]; 
depths= [procTbl.Raw_Depth_other_shark, procTbl.Second_Shark_Depth, procTbl.Third_Shark_Depth]; 

ns= length(Shark_Name); 
np= size(positions, 2); 

ci= ceil( linspace(1, size(cmap2, 1), ns+1) ); % define some colors

OS_code= ["AL", "TA", "TR", "YU"]; 
relative_positions= ["A", "B"];  % ahead or behind

% inspect to ensure key is standardized for sharks
% figure; histogram(categorical(other_sharks(:)))

% pre-allocate some big logical matrices
obs_and_other_behind= false([size(positions) ns]);
is_other_and_ahead=   false([size(positions) ns]);

obs_and_other_ahead=  false([size(positions) ns]);
is_other_and_behind=  false([size(positions) ns]);

% get broad lead-follow interactions from data
for s= 1:ns
    curr_shark= Shark == Shark_Name(s); 
    tmp_OS=  other_sharks == OS_code(s); 
    for p= 1:np
        obs_and_other_behind(:, p, s)= curr_shark & (positions(:, p) == "B" | strncmpi(positions(:, p), "B", 1)); 
        is_other_and_ahead(:, p, s)=   tmp_OS(:, p) & (positions(:, p) == "A" | strncmpi(positions(:, p), "A", 1));
        
        shark_leads(:, p)=  obs_and_other_behind(:, p, s) | is_other_and_ahead(:, p, s); 

        obs_and_other_ahead(:, p, s)=  curr_shark & (positions(:, p) == "A" | strncmpi(positions(:, p), "A", 1));
        is_other_and_behind(:, p, s)=  tmp_OS(:, p) & (positions(:, p) == "B" | strncmpi(positions(:, p), "B", 1));

        shark_follows(:, p)= obs_and_other_ahead(:, p, s) | is_other_and_behind(:, p, s); 
    end

    leads(s)= sum(shark_leads(:)); 
    follows(s)= sum(shark_follows(:));
end

tmp_matrix= [Shark Shark Shark]; 
follow_counts= cell(4, 1); 

lead_counts= follow_counts; 

% get counts and identity of leading shark's followers
for s= 1:ns
    tmp_dat1= other_sharks(squeeze(obs_and_other_behind(:, :, s)));
    tmp_dat2= tmp_matrix(squeeze(is_other_and_ahead(:, :, s))); 

    follow_counts{s}= categorical([tmp_dat1; tmp_dat2]);

    tmp_dat3= other_sharks(squeeze(obs_and_other_ahead(:, :, s)));
    tmp_dat4= tmp_matrix(squeeze(is_other_and_behind(:, :, s))); 

    lead_counts{s}= categorical([tmp_dat3; tmp_dat4]); 
end

[hc, cats]= cellfun(@histcounts, follow_counts, 'UniformOutput', false);
[hc2, cats2]= cellfun(@histcounts, lead_counts, 'UniformOutput', false);


% construct lead-follow matrix from data
final_matrix= zeros(4); 

for s= 1:ns
    this_shark= ismember(cats{s}, Shark_Name(s));
    ind= find(this_shark);
    ix_others= ~this_shark & ismember(cats{s}, Shark_Name); 
    ix= ~ismember(1:ns, s); 
    final_matrix(ix, s)= hc{s}(ix_others);

    leads(s)= leads(s) - hc{s}(ind);       % fix incorrect self-follows

    ts2= ismember(cats2{s}, Shark_Name(s));
    i2= find(ts2);
    follows(s)= follows(s) - hc2{s}(i2);   % fix incorrect self-leads
end

LF_matrix= [leads; follows]; 


%% Plot Figure 5 panels A and B

figure('Name', 'Interaction participation'); 
subplot(1, 2, 1);
barh((4:-1:1)', LF_matrix', 'stacked')
yticklabels(flipud(Shark_Name))
subplot(1, 2, 2);
circPercent(sum(LF_matrix, 1) ./ sum(LF_matrix(:)), 2, 'color', cmap2(ci(1:4), :), 'precision', 1);

LF_index= (LF_matrix(1, :) - LF_matrix(2, :)) ./ sum(LF_matrix, 1);

figure('Name', 'Lead-follow index'); hold on;
set(gcf, 'color', 'w')
yline(0)
scatter(1:ns, LF_index, 200, 1:ns, 'filled'); 
colormap(cmap2(ci(1:4), :))
title('lead-follow index')

% set box property to off and remove background color
ax= gca; 
set(ax, 'box', 'off', 'color', 'none', 'TickDir', 'out')
ax.XAxis.Visible= 'off';

xlim([0 ns+1])
ylabel('LFI')
grid on; 


%% Matrix supplied by Al

if useComputed == 1
    LF_obs_mat= final_matrix; 
else
    % recreation of Al's table
    LF_obs_mat= [0 363 262 94; 
                 386 0 186 90; 
                 338 213 0 68; 
                 148 87 99 0];
end

N= sum(LF_obs_mat(:));

% compare computed proportions with given
% computed= LF_matrix(1, :) ./ sum(LF_matrix(1, :)); 
% given= sum(LF_obs_mat, 1) ./ N; 

% figure; 
% circPercent([computed; given], 2, 'precision', 1);

% figure; 
% matrixTable(LF_obs_mat, Shark_Name, Shark_Name, 'XAxisLocation', 'top', ...
%             'dualX', string(sum(LF_obs_mat, 1)), 'dualY', string(sum(LF_obs_mat, 2)));
% figure; 
% matrixTable(original_data, Shark_Name, Shark_Name, 'XAxisLocation', 'top', ...
%             'dualX', string(sum(original_data, 1)), 'dualY', string(sum(original_data, 2)));


%% Chi-square 

row_totals= sum(LF_obs_mat, 2);
col_totals= sum(LF_obs_mat, 1); 
lf_totals=  sum(row_totals); 

lf_exp= ceil( row_totals * col_totals / lf_totals ); 

lf_stats= X2ind(LF_obs_mat, lf_exp, 0.05); 

LF_stat_mat= lf_stats.Computed;
zeroDiag= ~eye(size(LF_stat_mat));   % zero out diagonal
LF_stat_mat= LF_stat_mat .* zeroDiag;

overall_p= chi2cdf( sum(LF_stat_mat(:)), 9, 'upper' ); 

% p-values for individual leader-follower pairs
pchi= chi2cdf( LF_stat_mat, 9, 'upper' ); 

n_comps= numel(pchi) - sum(~zeroDiag(:)); 
p_thresh= 0.05 / n_comps; 

O_E= LF_obs_mat ./ lf_exp;   % observed:expected ratio 


%% Plot Figure 5 panels C and D

% make a digraph object
% the digraph object 'G' contains an Edges table (syntax --> 'G.Edges')
% that we populate below with the ratio of Observed:Expected, p-values, and
% indicators of statistical significance
G= digraph(LF_obs_mat, Shark_Name, 'omitselfloops'); 
G.Edges.LWidths= G.Edges.Weight / 80; 
G.Edges.RatioOE= O_E(zeroDiag); 
G.Edges.pVal= pchi(zeroDiag);
G.Edges.Sig= strings(size(pchi(zeroDiag))); 
G.Edges.Sig(pchi(zeroDiag) < p_thresh)= '***';
G.Edges.CData= -log(G.Edges.pVal);

H= flipedge(G);  % flip edges so that arrow direction indicates following

% plot matrix and graph network
figure('Name', 'LF pair comparisons'); 
subplot(1, 2, 1);
h= xcorrDotPlot(-log(pchi), 'full', 'dotscale', 150,  ...
                'rowlabels', Shark_Name, 'collabels', Shark_Name, ...
                'majorgrid', 0, 'minorgrid', 1, 'alpha', 1, ...
                'overlayvals', 'all', 'precision', 2);
colormap(cmap)

clim([-2.5 max(G.Edges.CData)]);


subplot(1, 2, 2);
colormap(cmap);
pG= plot(H, 'EdgeLabel', H.Edges.Sig);
pG.MarkerSize= 12;  
pG.NodeColor= 'k';  pG.NodeFontSize= 14;
pG.EdgeCData= H.Edges.CData;  pG.EdgeFontSize= 14;  pG.EdgeAlpha= 1;
pG.ArrowPosition= 0.9;  pG.ArrowSize= 12.5 * O_E(zeroDiag);
pG.LineWidth= H.Edges.LWidths;
set(gcf,'color','w');
clim([-2.5 max(G.Edges.CData)])

axis equal