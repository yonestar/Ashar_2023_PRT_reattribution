basedir = '/Users/yoni/Repositories/Ashar_2023_PRT_reattribution/';
figdir = fullfile(basedir, 'figures');

% Note on providence: reidentify_data.m takes IPQlastitem_final_deid categories reconciled.csv, generates IPQlastitem categories.csv
% merge_with_metadata.m takes IPQlastitem categories.csv, generates 'IPQlastitem_final MERGED LONG CATEGORIES.csv'
% Which is now ready for analyses

cats_long = readtable(fullfile(basedir, 'data', 'manual_codings', 'IPQlastitem_final MERGED LONG CATEGORIES.csv'));
cats_wide = readtable(fullfile(basedir, 'data', 'manual_codings', 'IPQlastitem_final MERGED WIDE CATEGORIES.csv')); % for mediation

cats_long.time = strcmp(cats_long.redcap_event_name, 't2_arm_1') + 1;

head(cats_long)

%% confirm unique categories
unique(cats_long.cat1)'
unique(cats_long.cat2)'
unique(cats_long.cat3)'

%% combine low prevalence categories

% put job into activity
% put diet into other
% put weight into physio
% put perinatal into activity
% combine posture with physio

for i=6:8 %cols for the three cats
    
    wh = strcmp(cats_long{:,i}, 'job');
    if (sum(wh)>0), cats_long(wh,i) = repmat({'activity'}, sum(wh), 1); end
    
    wh = strcmp(cats_long{:,i}, 'diet');
    if (sum(wh)>0), cats_long(wh,i) = repmat({'other'}, sum(wh), 1); end
    
    wh = strcmp(cats_long{:,i}, 'weight');
    if (sum(wh)>0), cats_long(wh,i) = repmat({'physiology (non-spinal)'}, sum(wh), 1); end

    wh = strcmp(cats_long{:,i}, 'posture');
    if (sum(wh)>0), cats_long(wh,i) = repmat({'physiology (non-spinal)'}, sum(wh), 1); end

    wh = strcmp(cats_long{:,i}, 'perinatal');
    if (sum(wh)>0), cats_long(wh,i) = repmat({'activity'}, sum(wh), 1); end
    
    
end

%% convert to categorical
cats_long.cat1 = categorical(cats_long.cat1);
cats_long.cat2 = categorical(cats_long.cat2);
cats_long.cat3 = categorical(cats_long.cat3);

cats_wide.cat1_t1_arm_1 = categorical(cats_wide.cat1_t1_arm_1);
cats_wide.cat2_t1_arm_1 = categorical(cats_wide.cat2_t1_arm_1);
cats_wide.cat3_t1_arm_1 = categorical(cats_wide.cat3_t1_arm_1);

cats_wide.cat1_t2_arm_1 = categorical(cats_wide.cat1_t2_arm_1);
cats_wide.cat2_t2_arm_1 = categorical(cats_wide.cat2_t2_arm_1);
cats_wide.cat3_t2_arm_1 = categorical(cats_wide.cat3_t2_arm_1);

%% plot treemaps 

mycats = cats_long(cats_long.group<5 & cats_long.time==1,:);
plot_treemap([mycats.cat1; mycats.cat2; mycats.cat3], 'treemap T1 PRT')  

mycats = cats_long(cats_long.group==1 & cats_long.time==2,:);
plot_treemap([mycats.cat1; mycats.cat2; mycats.cat3], 'treemap T2 PRT') 

% combined controls
mycats = cats_long(cats_long.group>1 & cats_long.time==1,:); plot_treemap([mycats.cat1;
mycats.cat2; mycats.cat3], 'treemap T1 controls');

mycats = cats_long(cats_long.group>1 & cats_long.time==2,:); plot_treemap([mycats.cat1;
mycats.cat2; mycats.cat3], 'treemap T2 controls')

% Placebo and UC shown separately
% mycats = cats(cats.group==2 & cats.time==1,:);
% plot_treemap([mycats.cat1; mycats.cat2; mycats.cat3], 'treemap T1 Placebo')  
% 
% mycats = cats(cats.group==2 & cats.time==2,:);
% plot_treemap([mycats.cat1; mycats.cat2; mycats.cat3], 'treemap T2 Placebo') 
% 
% mycats = cats(cats.group==3 & cats.time==1,:);
% plot_treemap([mycats.cat1; mycats.cat2; mycats.cat3], 'treemap T1 UC')  
% 
% mycats = cats(cats.group==3 & cats.time==2,:);
% plot_treemap([mycats.cat1; mycats.cat2; mycats.cat3], 'treemap T2 UC') 


%% make a stacked bar graph

% build the percentages data structure
% percentages of attribution categoriesone row for each group at each timepoint (each row) 
[percentages(6,:), rawcounts(6,:)] = get_cat_percentages(cats_long, cats_long.group==1 & cats_long.time==1);
[percentages(5,:), rawcounts(5,:)] = get_cat_percentages(cats_long, cats_long.group==1 & cats_long.time==2);
[percentages(4,:), rawcounts(4,:)] = get_cat_percentages(cats_long, cats_long.group==2 & cats_long.time==1);
[percentages(3,:), rawcounts(3,:)] = get_cat_percentages(cats_long, cats_long.group==2 & cats_long.time==2);
[percentages(2,:), rawcounts(2,:)] = get_cat_percentages(cats_long, cats_long.group==3 & cats_long.time==1);
[percentages(1,:), rawcounts(1,:)] = get_cat_percentages(cats_long, cats_long.group==3 & cats_long.time==2);

labels = categories([cats_long.cat1; cats_long.cat2; cats_long.cat3]);
labels{11} = 'spinal'; labels{8} = 'physio'; labels{4} = 'genetic'; labels{9} = 'psych'; labels{6} = 'neglect';    

% make brain, psych, and stress the first 3
% swap brain with activities
[percentages, labels, rawcounts] = swapsies(percentages, labels, rawcounts, 2, 3);
[percentages, labels, rawcounts] = swapsies(percentages, labels, rawcounts, 1, 9);
[percentages, labels, rawcounts] = swapsies(percentages, labels, rawcounts, 3, 12);
[percentages, labels, rawcounts] = swapsies(percentages, labels, rawcounts, 5, 11);
[percentages, labels, rawcounts] = swapsies(percentages, labels, rawcounts, 9, 12);

create_figure('stackedbar'); 
h = barh([1 2 4 5 7 8], percentages*100, 'stacked');
set(gca, 'FontSize', 18)
xlim([0 130]) % to make room for legend
legend(labels, 'FontSize', 18, 'Location', 'east')
xlabel('% attributions')

colors = brewermap(numel(labels), 'Paired'); % qualitative

% switch colors for beauty
tmp = colors(3,:);
colors(3,:) = colors(4,:);
colors(4,:) = tmp;

tmp = colors(1,:);
colors(1,:) = colors(2,:);
colors(2,:) = tmp;

tmp = colors(11,:);
colors(11,:) = colors(10,:);
colors(10,:) = tmp;

% Set colors for each bar
for i = 1:numel(labels)
    h(i).FaceColor = colors(i,:);
end

% label
yticks = 1:8;
yticklabels = {'PRT Pre', 'PRT Post', '', 'Placebo Pre', 'Placebo Post', '', 'Usual Care Pre', 'Usual Care Post'};  % Replace with your desired labels
set(gca, 'YTick', yticks, 'YTickLabel', flip(yticklabels), 'XTick', 0:20:100);

%% export stacked bar in vector format
% adjust size label and legend position manually then execute below

% Save the figure using the 'print' function
delete(fullfile(figdir,'stacked_bar.pdf'));
print(gcf, fullfile(figdir,'stacked_bar.pdf'), ['-d', 'pdf'], '-r300', '-bestfit');

%% make supp table with values
clc
fliplr(round(percentages'*100))
fliplr(rawcounts')
yticklabels

%% make a stacked bar graph -- COMBINED CONTROLS

% combine final rows of above
percentages_combined(3:4,:) = percentages(5:6,:);
percentages_combined(1,:) = mean(percentages([1 3],:));
percentages_combined(2,:) = mean(percentages([2 4],:));

create_figure('stackedbar 2'); 
h = barh([1 2 4 5], percentages_combined*100, 'stacked');
set(gca, 'FontSize', 18)
xlim([0 130]) % to make room for legend
legend(labels, 'FontSize', 18, 'Location', 'east')
xlabel('% attributions')

colors = brewermap(numel(labels), 'Paired'); % qualitative

% switch colors for beauty
tmp = colors(3,:);
colors(3,:) = colors(4,:);
colors(4,:) = tmp;

tmp = colors(1,:);
colors(1,:) = colors(2,:);
colors(2,:) = tmp;

tmp = colors(11,:);
colors(11,:) = colors(10,:);
colors(10,:) = tmp;

% Set colors for each bar
for i = 1:numel(labels)
    h(i).FaceColor = colors(i,:);
end

% label
yticks = 1:5;
yticklabels = {'PRT Pre', 'PRT Post', '', 'Controls Pre', 'Controls Post'};  % Replace with your desired labels
set(gca, 'YTick', yticks, 'YTickLabel', flip(yticklabels), 'XTick', 0:20:100);


%% create MBA scores for  psych/brain/stress and write back to file

% make a score for # of psych/brain
cats_long.mindbrain_attr_score_from_cats = cat_is_mindbrain(cats_long.cat1) + cat_is_mindbrain(cats_long.cat2) + cat_is_mindbrain(cats_long.cat3);
writetable(cats_long, fullfile(basedir, 'data', 'manual_codings', 'IPQlastitem_final MERGED LONG CATEGORIES.csv'));

% import score into cats_wide -- verified correct
cats_wide.mindbrain_attr_score_from_cats_t1_arm_1 = cat_is_mindbrain(cats_wide.cat1_t1_arm_1) + cat_is_mindbrain(cats_wide.cat2_t1_arm_1) + cat_is_mindbrain(cats_wide.cat3_t1_arm_1);
cats_wide.mindbrain_attr_score_from_cats_t2_arm_1 = cat_is_mindbrain(cats_wide.cat1_t2_arm_1) + cat_is_mindbrain(cats_wide.cat2_t2_arm_1) + cat_is_mindbrain(cats_wide.cat3_t2_arm_1);

% dont write it out, b/c <undefined> becomes '<undefined>' string which
% causes errors when attempting to read it back in later
% writetable(cats_wide, fullfile(basedir, 'data', 'manual_codings', 'IPQlastitem_final MERGED WIDE CATEGORIES.csv'));

%% percent mind/brain at pre and post
clc
wh = cats_long.time == 2 & cats_long.group==1;
sum(cats_long.mindbrain_attr_score_from_cats(wh))
sum(cats_long.mindbrain_attr_score_from_cats(wh)) / (numel(cats_long.mindbrain_attr_score_from_cats(wh)) * 3)

wh = cats_long.time == 2 & cats_long.group==2;
sum(cats_long.mindbrain_attr_score_from_cats(wh))
sum(cats_long.mindbrain_attr_score_from_cats(wh)) / (numel(cats_long.mindbrain_attr_score_from_cats(wh)) * 3)

wh = cats_long.time == 2 & cats_long.group==3;
sum(cats_long.mindbrain_attr_score_from_cats(wh))
sum(cats_long.mindbrain_attr_score_from_cats(wh)) / (numel(cats_long.mindbrain_attr_score_from_cats(wh)) * 3)


%% plot group by time, MBA scores

[h1,h2,h3] = plot_olp4cbp_group_by_time(cats_long, 'mindbrain_attr_score_from_cats');

% save figures
close(h3)
figure(h1)
set(h1, 'Position', [ 284   247   446   313])

ylabel({'# of ''psychological'', ''stress'',', 'or ''brain'' attributions (0-3)'})
legend('Location', 'NW')
print_pdf(h1, fullfile(figdir, 'attributions_grp_by_time_from_categories.pdf'))
set(gca, 'FontSize', 20);

figure(h2)
set(h2, 'Position', [ 784    182   464   378])
set(gca, 'FontSize', 20)
xlabel('Δ Pain, Pre-to-post-tx');
ylabel({'# ''psych'', ''stress'', or ''brain''', 'Δ Pre-to-post-tx'})
%ylabel({'Mind-brain attribution', 'Δ Pre-to-post-tx'});
%legend({'PRT' 'Placebo' 'No Tx'});
legend off
print_pdf(h2, fullfile(figdir, 'attributions_ind_diffs_from_categories.pdf'))

%% Effect of group on change in mind brain
clc
cats_delta = create_metadata_delta(cats_long);
cats_delta.pain_avg_delta_z = zscore(cats_delta.pain_avg_delta);
cats_delta.mindbrain_attr_score_from_cats_delta_z = zscore(cats_delta.mindbrain_attr_score_from_cats_delta);

fitglm(cats_delta, 'mindbrain_attr_score_from_cats_delta_z ~ group + age + gender', 'CategoricalVars', 'group')

%% non-param effect size estimates for PRT vs control on MBA scores -- hedges g
cats_delta = create_metadata_delta(cats_long);
delta_attr = cats_delta.mindbrain_attr_score_from_cats_delta;
clc

stats = mes(delta_attr(cats_delta.group==1), delta_attr(cats_delta.group==2), 'hedgesg')
stats = mes(delta_attr(cats_delta.group==1), delta_attr(cats_delta.group==3), 'hedgesg')
stats = mes(delta_attr(cats_delta.group==3), delta_attr(cats_delta.group==2), 'hedgesg')

%% Assocation of pre to post change in MBA scores and change in pain intensity  
clc
cats_delta = create_metadata_delta(cats_long);

cats_delta.pain_avg_delta_z = zscore(cats_delta.pain_avg_delta);
cats_delta.mindbrain_attr_score_from_cats_delta_z = zscore(cats_delta.mindbrain_attr_score_from_cats_delta);

fitglm(cats_delta, 'pain_avg_delta_z ~ mindbrain_attr_score_from_cats_delta_z * group + age + gender', 'CategoricalVars', 'group')

% follow-up: correlation of simple change scores: change in MBAS and pain within group or whole sample
clc, 
wh = cats_delta.group>0; % full sample
wh = cats_delta.group==1; % just PRT

% delta_attr created in cell above
[r,p] = corr(cats_delta.pain_avg_delta(wh), delta_attr(wh))
[r,p] = corr(cats_delta.tsk11_delta(wh), delta_attr(wh))
[r,p] = corr(cats_delta.pcs_delta(wh), delta_attr(wh))

%% Assocation of pre post change in MBA scores and change in pain intensity to 1 yr
clc

cats_wide.pain_1yr_delta_z = scale(cats_wide.pain_avg_x12_month_follow_up_arm_1 - cats_wide.pain_avg_baseline);
cats_wide.mindbrainscore_delta_z = scale(cats_wide.mindbrain_attr_score_from_cats_t2_arm_1 - cats_wide.mindbrain_attr_score_from_cats_t1_arm_1);

fitglm(cats_wide, 'pain_1yr_delta_z ~ mindbrainscore_delta_z * group + age + gender', 'CategoricalVars', 'group')

% follow-up: correlation of simple change scores: change in MBAS and pain within group or whole sample
wh = cats_delta.group>0;
wh = cats_delta.group==1;

[r,p] = corr(cats_wide.pain_1yr_delta_z(wh), cats_wide.mindbrainscore_delta_z(wh), 'rows', 'complete')

%% GLM of change scores and TSK11, standardized 
clc
cats_delta.tsk11_delta_z = zscore(cats_delta.tsk11_delta);
fitglm(cats_delta, 'tsk11_delta_z ~ mindbrain_attr_score_from_cats_delta_z * group + age + gender', 'CategoricalVars', 'group')

%% GLM of change scores and SOPA, standardized
clc
cats_delta.sopa_emo_delta_z = zscore(cats_delta.sopa_emo_delta);
fitglm(cats_delta, 'sopa_emo_delta_z ~ mindbrain_attr_score_from_cats_delta_z * group + age + gender', 'CategoricalVars', 'group')

%% GLM of change scores and PCS, standardized 
clc
cats_delta.pcs_delta_z = zscore(cats_delta.pcs_delta);
fitglm(cats_delta, 'pcs_delta_z ~ mindbrain_attr_score_from_cats_delta_z * group + age + gender', 'CategoricalVars', 'group')

%% distribution of baseline mb scores
x = cats_long.mindbrain_attr_score_from_cats(cats_long.time==1);
mean(x), median(x), min(x), max(x)


%% assocation between gender and baseline attribution
% perm test since non-normal
wh = cats_long.time==1;
whg = cats_long.gender==1;

figure; 
histogram(cats_long.mindbrain_attr_score_from_cats(wh & whg))
hold on,
histogram(cats_long.mindbrain_attr_score_from_cats(wh & ~whg))

[p, observeddifference, effectsize] = permutationTest(cats_long.mindbrain_attr_score_from_cats(wh & whg), cats_long.mindbrain_attr_score_from_cats(wh & ~whg), 10000, 'plotresult', 0)

ttest2(cats_long.mindbrain_attr_score_from_cats(wh & whg), cats_long.mindbrain_attr_score_from_cats(wh & ~whg))

%% assocation between baseline attribution and X
wh = cats_long.time==1;
Y = cats_long.mindbrain_attr_score_from_cats(wh);


X = cats_long.backpain_length(wh);
%X = cats.age(wh);
X = cats_long.pain_avg(wh);
% 
% X = cats_long.tsk11(wh);
% X = cats_long.pcs(wh);
% 
X = cats_long.sopa_emo(wh);

clc
if 0
    figure, scatterhist(X, Y);
    set(gca, 'FontSize', 24), ylabel('Baseline attribution');xlabel('X'); lsline
end
[B,STATS] = robustfit(zscore(Y), zscore(X));
B(2), STATS.p(2)

[r,p]=corr(X, Y, 'Type', 'Spearman')


%% mediation: X is PRT vs control, M is MBA scores at mid-tx, Y is pain intensity at 1 year
clc

% rename a variable, as expected by the sub-function
cats_wide.mindbrain_attr_score_from_cats_baseline = cats_wide.mindbrain_attr_score_from_cats_t1_arm_1;

mediation_to_follow_up_combined_controls(cats_wide, 'mindbrain_attr_score_from_cats', 'pain_avg')

%% subfunctions

% for swapping orders in stacked bar graph
function [percentages, labels, rawcounts] = swapsies(percentages, labels, rawcounts, ind1, ind2)

    tmp = percentages(:,ind1);
    percentages(:,ind1) = percentages(:,ind2);
    percentages(:,ind2) = tmp;
    
    tmp = labels(ind1);
    labels(ind1) = labels(ind2);
    labels(ind2) = tmp;

    tmp = rawcounts(:,ind1);
    rawcounts(:,ind1) = rawcounts(:,ind2);
    rawcounts(:,ind2) = tmp;
end

function [percentages, rawcount] = get_cat_percentages(cats, wh)
    mycats = cats(wh,:);
    dat = [mycats.cat1; mycats.cat2; mycats.cat3];
    percentages = countcats(dat) / numel(dat);
    rawcount = countcats(dat);
end

% input: array of categories
% output: 0/1 whether each category counts as mind/brain
function is_mb = cat_is_mindbrain(cat_array)
    
    is_mb = (cat_array == 'psychological') + (cat_array == 'brain') + (cat_array == 'stress');
    is_mb(isundefined(cat_array)) = NaN;
end   

function plot_pie(dat, fig_name)
    basedir = '/Users/yoni/Repositories/Ashar_2023_PRT_reattribution/';
    figdir = fullfile(basedir, 'figures');

    figure
    h = pie(dat);
    set(gca, 'visible', 'off')
    set(h(2:2:end),'FontSize',20);

    set(gcf, 'Position', [   693   554   624   319], 'Color', 'white');
    
    % printing doesn't work well. 1) manually adjust labels so they look
    % good, 2) open in print preview, 3) change to 'auto' and save to PDF
    
%    saveas(gcf, fullfile(figdir, [fig_name '.pdf']));
    %print(gcf, '-dpdf', fullfile(figdir, [fig_name '.pdf']), '-fillpage');

end

function plot_treemap(dat, fig_name)
    basedir = '/Users/yoni/Repositories/Ashar_2023_PRT_reattribution/';
    figdir = fullfile(basedir, 'figures');

    figure
    colors = brewermap(numel(categories(dat)), 'set3'); % qualitative

    counts = countcats(dat);
    
    % switch colors 6 and 9, 8 and 10, for beauty
    tmp = colors(6,:);
    colors(6,:) = colors(9,:);
    colors(9,:) = tmp;
    
    tmp = colors(8,:);
    colors(8,:) = colors(10,:);
    colors(10,:) = tmp;
    
    labels = categories(dat);
    labels{11} = 'spinal'; labels{8} = 'physio'; labels{4} = 'genetic'; labels{9} = 'psych'; labels{6} = 'neglect';
    
    % remove any empty category
    if any(counts==0)
        wh = counts==0;
        counts(wh) = [];
        labels(wh) = [];
        colors(wh,:) = [];
    end
    
    rects = treemap(counts);

    for i=1:length(labels)
        labels{i} = sprintf('%s (%i%%)', labels{i}, round(100*counts(i)/sum(counts)));
    end
    
    plotRectangles(rects,labels, colors);
    outline(rects)
    axis([-0.01 1.01 -0.01 1.01])

    
    % saveas(gcf, fullfile(figdir, [fig_name '.pdf']));
    % NO GOOD print(gcf, '-dpdf', fullfile(figdir, [fig_name '.pdf']), '-fillpage');

end

% mediation for M at post-tx, Outcome at all timepoints, using change
% scores
function STATS =  mediation_to_follow_up_combined_controls(outcomes_wide, Mname, Yname)

    outcome = get_outcome_data_from_final_acute_outcomes_wide_prt(outcomes_wide, Yname, 'noplot', 'nogpml', '1yearnoweekly'); % substracts baseline by default
    mymediator = get_outcome_data_from_final_acute_outcomes_wide_prt(outcomes_wide, Mname, 'noplot', 'nogpml', 'pre_post_only'); % substracts baseline by default
    posttx_ind = find(strcmp(outcome.xlab, ' Post-tx'));
    
    X = [ones(size(outcome.datmat_prt, 1), 1); zeros(size(outcome.datmat_control, 1), 1); zeros(size(outcome.datmat_pla, 1), 1)];
    outcome_datmat = [outcome.datmat_prt; outcome.datmat_control; outcome.datmat_pla];
    mediator_datmat = [mymediator.datmat_prt; mymediator.datmat_control; mymediator.datmat_pla];

    M = mediator_datmat(:, posttx_ind); % mediator at post-tx
    cov = [outcomes_wide.age outcomes_wide.gender]; % covariates
    
    clc, STATS = {};
    % loop through outcome timepoints (Y). X and M are same for all
    % iteration
    for i=6 % from post-tx to last timepiont

        print_header(['M = Delta ' Mname ' from baseline to post-tx'], ['Y = ' Yname ' at ' outcome.xlab{i + posttx_ind - 1}])
        Y = outcome_datmat(:, i + posttx_ind - 1); % outcome at given timepoint

        [~, stats] = mediation(X,scale(Y),scale(M), 'verbose',  'boot', 'bootsamples', 10000, 'covs', cov, 'doCIs');
                 
        STATS.a(i, 1) = stats.z(1); % mean(1) ./ stats.ste(1);
        STATS.a_p(i, 1) = stats.p(1);

        STATS.b(i, 1) = stats.z(2); %./ stats.ste(2);
        STATS.b_p(i, 1) = stats.p(2);

        STATS.ab(i, 1) = stats.z(5); %./ stats.ste(5);
        STATS.ab_p(i, 1) = stats.p(5);

        STATS.perc_mediated(i, 1) = 100 .* (stats.mean(4) - stats.mean(3)) ./ stats.mean(4);
    end
end
