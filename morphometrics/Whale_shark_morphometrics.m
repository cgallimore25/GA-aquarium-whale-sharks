%% Analyze whale shark sizes and weights

% Key
% PCL - Pre-caudal length
% PCP - Pre-caudal pit length
% TL -  Total length
% FL -  Fork length 
% PD1 - Pre-first dorsal length; snout tip to first dorsal fin origin; 
%       most appropriate estimate of TL
% HW -  Head width (Eye to eye)
% MOW - Mouth width

% others from literature
% PD2 - Pre-second dorsal length
% PP2 - Pre-pelvic length
% PAL - Pre-anal length
% CDM - Length of dorsal caudal margin
% INO - Interorbital space

% D1A - First Dorsal Fin anterior margin (top/leading edge)
% D1B - First Dorsal-Fin Base Length
% D1H - First Dorsal-Fin Height

% see Elasmobranch husbandry manual II, Ch. 2, pg 18, Table 2,
% Rhincodon typus (Smith, 1828) &
% Rogers et al., 2017, Marine Biology

% correlation between log(PD1) and log(TL) and Body Mass (BM) estimation
% equation (pg 19)


%% Parameters

plotALLmorphs= 0; 

% make a colormap
cmap= cmocean('ice');   

% specify file
fileName= 'WS_Morphometrics_all.xlsx'; 

opts= detectImportOptions( fileName );

opts= setvartype(opts, {'Shark', 'Date', 'MeasurementType', 'Unit'}, 'string');      

dataTbl = readtable(fileName, opts);  

dates= dataTbl.Date;  
dates_dt= datetime(dates, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss', 'TimeZone', 'America/New_York');

Shark= upper( dataTbl.Shark );

units= dataTbl.Unit; 

meas_type= dataTbl.MeasurementType; 
quant= dataTbl.Quantity; 

u_mt= unique(meas_type); 

u_sharks= unique(Shark);
ns= length(u_sharks); 

des_mt= ["D1A", "PD1", "FL", "TL", "PCL", "HW", "MOW"]; 
nm= length(des_mt); 

yrs= year(dates_dt); 
u_y= unique(yrs); 
des_y= u_y(1:4); 
ny= length(des_y);

morph_mat= zeros(ns, ny, nm); 
morph_std= nan(size(morph_mat)); 
n_meas= morph_mat; 

n_el= ns * ny * nm; 

c_idx= linspace(1, size(cmap, 1), ns); % color indices
cols= cmap(c_idx, :);                  % define some colors


c= 1; 

for s= 1:ns
    sv= Shark == u_sharks(s);
    for y= 1:ny
        yv= yrs == des_y(y);
        % get weight
        wv= meas_type == "Weight"; 
        w= sv & yv & wv; 
        if sum(w) == 0
            weight_mat(s, y)= NaN;
        elseif sum(w) > 1
            weight_mat(s, y)= mean(quant(w));
        else
            weight_mat(s, y)= quant(w);
        end
        % get length metrics
        for m= 1:nm
            mv= meas_type == des_mt(m); 
            val= sv & yv & mv; 

            if sum(val) == 0
                morph_mat(s, y, m)= NaN; 
            elseif sum(val) > 1
                morph_mat(s, y, m)= mean(quant(val)); 
                morph_std(s, y, m)= std(quant(val)); 
            else
                morph_mat(s, y, m)= quant(val); 
            end

            sum_vec(c)= sum(val); 

            n_meas(s, y, m)= sum(val); 

            c= c + 1; 
        end
    end
end


if plotALLmorphs  % plot for all measures

    for m= 1:nm
        dat= squeeze(morph_mat(:, :, m)); 
        figure; plot(dat', '--o'); 
        title(des_mt(m))
        xlim([0 5])
        xticks(0:5)
        xticklabels([""; des_y; ""]); 
        xlabel('year')
        ylabel('cm')
        legend(u_sharks)
    end
end


% Process by averaging across the years 2008 and 2009 (length), and take 
% average weights for the year 2009

% Make either a summarizing radial plot or bivariate plot of TL by weight

TL=  squeeze(morph_mat(:, :, 4)); 
PD1= squeeze(morph_mat(:, :, 2)); 
FL=  squeeze(morph_mat(:, :, 3)); 
PCL= squeeze(morph_mat(:, :, 5)); 

av_TL= sum([TL(:, 2) TL(:, 3)], 2) ./ 2;
av_FL= sum([FL(:, 2) FL(:, 3)], 2) ./ 2;
av_PD1= sum([PD1(:, 2) PD1(:, 3)], 2) ./ 2;
av_PCL= sum([PCL(:, 2) PCL(:, 3)], 2) ./ 2;

TL0= TL(:, 1)/100; 
lengths_ft= convlength(sum([TL(:, 2) TL(:, 3)], 2) ./ (100 * 2), 'm', 'ft' ); 
lengths_photog_ft= [26.7; 24.0; 31.2; 26.6]; 
lengths_photog_m= convlength(lengths_photog_ft, 'ft', 'm' ); 


% plot relationship of Pre-first dorsal-fin length (PD1) & Total length (TL)
% from Figure 2, Elasmobranch husbandry manual II, Ch. 2, pg 19
figure; 
scatter(log10(av_PD1), log10(av_TL))
xlim([2.15 2.6])
ylim([2.5 3])



%% plot individual concordant metrics

% colors and line/marker styles
c2= distinguishable_colors(8);
c2= c2(5:end, :); 
darkF= [1.85; 1.2; 1]; % darkness factor

linestyles= ["-", "--", "-."];
mrkrstyles= ["o", "square", "^"];

concordant_morphs= cat(3, TL, FL, PCL); 
m_tmp= permute(concordant_morphs, [2 3 1]);

yL= [350 700; 300 600; 260 510]; 
titls= ["total length", "fork length", "pre-caudal length", "combined"]; 

figure; 
for m= 1:3
    subplot(1, 4, m);
    tc= c2 ./ darkF(m);
    lspec= strcat(linestyles(m), mrkrstyles(m)); 
    for s= 1:ns
        sh_m= squeeze(m_tmp(:, :, s));
        plot(sh_m(:, m), lspec, 'Color', tc(s, :)); hold on
    end
    if m == 1
        ylabel('cm')
        xlabel('year')
        legend(u_sharks)
    end
    title(titls(m)); grid on
    xlim([0 5])
    xticks(0:5)
    xticklabels([""; des_y; ""]); 
    ylim(yL(m, :))
end


% plot max normalized (0-1) concordant morphometrics
mn_concordant_morphs= cat(3, max_normall(TL), max_normall(FL), max_normall(PCL)); 
mn_tmp= permute(mn_concordant_morphs, [2 3 1]); % put shark on page dim

for s= 1:ns
    for m= 1:3
        subplot(1, 4, 4)
        tc= repmat(c2(s, :), 3, 1) ./ darkF; 
        sh_mn= squeeze(mn_tmp(:, :, s)); 
        lspec= strcat(linestyles(m), mrkrstyles(m)); 
        plot(sh_mn(:, m), lspec, 'Color', tc(m, :)); hold on
    end
end
title(titls(4)); grid on
ylabel('normalized metric')
xlim([0 5])
xticks(0:5)
xticklabels([""; des_y; ""]); 
ylim([-.15 1.15])


%% Spider plot

c_idx2= ceil( linspace(1, size(cmap, 1), ns+1) ); % define some colors

morphometrics= [av_TL, av_FL, av_PCL, av_PD1, weight_mat(:, 4)];

spider_ax= [470, 388, 335, 198, 1328; 626, 537, 455, 258, 2000];

figure('Name', 'Shark morphmetrics');    clear sp2
sp2= spider_plot_class(morphometrics); 

sp2.AxesPrecision = 0;
sp2.AxesLabels= {'TL', 'FL', 'PCL', 'PD1', 'Weight'}; 
sp2.AxesLimits = spider_ax;
sp2.Color= cmap(c_idx2(1:4), :);
sp2.LegendLabels = cellstr(u_sharks)';
sp2.FillOption = 'on';
sp2.FillTransparency =  0.2;  % change to 1 for figure export


%% Scale bars for Figure 1 a-b cartoon schematic

xc= [1:(ny+1); 1:(ny+1)];
yc= [zeros(1, ny+1); av_TL' 100];

figure; 
plot(xc, yc, '-k', 'LineWidth', 2)
xlim([0 ny+2])
