basedir = '/Users/yoni/Repositories/Ashar_2023_PRT_reattribution/';
figdir = fullfile(basedir, 'figures');

% Note on providence: reidentify_data.m takes IPQlastitem_final_deid categories reconciled.csv, generates IPQlastitem categories.csv
% merge_with_metadata.m takes IPQlastitem categories.csv, generates 'IPQlastitem_final MERGED LONG CATEGORIES.csv'
% Which is now ready for analyses

cats = readtable(fullfile(basedir, 'data', 'manual_codings', 'IPQlastitem_final MERGED LONG CATEGORIES.csv'));
cats_wide = readtable(fullfile(basedir, 'data', 'manual_codings', 'IPQlastitem_final MERGED WIDE CATEGORIES.csv')); % for mediation

cats.time = strcmp(cats.redcap_event_name, 't2_arm_1') + 1;

head(cats)

%% confirm unique categories

unique(cats.cat1)'
unique(cats.cat2)'
unique(cats.cat3)'

%% combine low prevalence categories

% put job into activity
% put diet into other
% put weight into physio
% put perinatal into activity
% combine posture with physio

for i=6:8 %cols for the three cats
    
    wh = strcmp(cats{:,i}, 'job');
    if (sum(wh)>0), cats(wh,i) = repmat({'activity'}, sum(wh), 1); end
    
    wh = strcmp(cats{:,i}, 'diet');
    if (sum(wh)>0), cats(wh,i) = repmat({'other'}, sum(wh), 1); end
    
    wh = strcmp(cats{:,i}, 'weight');
    if (sum(wh)>0), cats(wh,i) = repmat({'physiology (non-spinal)'}, sum(wh), 1); end

    wh = strcmp(cats{:,i}, 'posture');
    if (sum(wh)>0), cats(wh,i) = repmat({'physiology (non-spinal)'}, sum(wh), 1); end

    wh = strcmp(cats{:,i}, 'perinatal');
    if (sum(wh)>0), cats(wh,i) = repmat({'activity'}, sum(wh), 1); end
    
    
end

%% convert to categorical
cats.cat1 = categorical(cats.cat1);
cats.cat2 = categorical(cats.cat2);
cats.cat3 = categorical(cats.cat3);

cats_wide.cat1_t1_arm_1 = categorical(cats_wide.cat1_t1_arm_1);
cats_wide.cat2_t1_arm_1 = categorical(cats_wide.cat2_t1_arm_1);
cats_wide.cat3_t1_arm_1 = categorical(cats_wide.cat3_t1_arm_1);

cats_wide.cat1_t2_arm_1 = categorical(cats_wide.cat1_t2_arm_1);
cats_wide.cat2_t2_arm_1 = categorical(cats_wide.cat2_t2_arm_1);
cats_wide.cat3_t2_arm_1 = categorical(cats_wide.cat3_t2_arm_1);

%% plot treemaps 

mycats = cats(cats.group<5 & cats.time==1,:);
plot_treemap([mycats.cat1; mycats.cat2; mycats.cat3], 'treemap T1 PRT')  

mycats = cats(cats.group==1 & cats.time==2,:);
plot_treemap([mycats.cat1; mycats.cat2; mycats.cat3], 'treemap T2 PRT') 

% combined controls
mycats = cats(cats.group>1 & cats.time==1,:); plot_treemap([mycats.cat1;
mycats.cat2; mycats.cat3], 'treemap T1 controls');

mycats = cats(cats.group>1 & cats.time==2,:); plot_treemap([mycats.cat1;
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

%% create MBA scores for  psych/brain/stress

% make a score for # of psych/brain
cats.mindbrain_attr_score_from_cats = cat_is_mindbrain(cats.cat1) + cat_is_mindbrain(cats.cat2) + cat_is_mindbrain(cats.cat3);
writetable(cats, fullfile(basedir, 'manual_codings', 'IPQlastitem_final MERGED LONG CATEGORIES.csv'));

% import score into cats_wide -- verified correct
cats_wide.mindbrain_attr_score_from_cats_t1_arm_1 = cat_is_mindbrain(cats_wide.cat1_t1_arm_1) + cat_is_mindbrain(cats_wide.cat2_t1_arm_1) + cat_is_mindbrain(cats_wide.cat3_t1_arm_1);
cats_wide.mindbrain_attr_score_from_cats_t2_arm_1 = cat_is_mindbrain(cats_wide.cat1_t2_arm_1) + cat_is_mindbrain(cats_wide.cat2_t2_arm_1) + cat_is_mindbrain(cats_wide.cat3_t2_arm_1);

%% percent mind/brain at pre and post
clc
wh = cats.time == 2 & cats.group==1;
sum(cats.mindbrain_attr_score_from_cats(wh))
sum(cats.mindbrain_attr_score_from_cats(wh)) / (numel(cats.mindbrain_attr_score_from_cats(wh)) * 3)

wh = cats.time == 2 & cats.group==2;
sum(cats.mindbrain_attr_score_from_cats(wh))
sum(cats.mindbrain_attr_score_from_cats(wh)) / (numel(cats.mindbrain_attr_score_from_cats(wh)) * 3)

wh = cats.time == 2 & cats.group==3;
sum(cats.mindbrain_attr_score_from_cats(wh))
sum(cats.mindbrain_attr_score_from_cats(wh)) / (numel(cats.mindbrain_attr_score_from_cats(wh)) * 3)


%% plot group by time

[h1,h2,h3] = plot_olp4cbp_group_by_time(cats, 'mindbrain_attr_score_from_cats');

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
cats_delta.pain_avg_delta_z = zscore(cats_delta.pain_avg_delta);
cats_delta.mindbrain_attr_score_from_cats_delta_z = zscore(cats_delta.mindbrain_attr_score_from_cats_delta);

fitglm(cats_delta, 'mindbrain_attr_score_from_cats_delta_z ~ group + age + gender', 'CategoricalVars', 'group')

%% non-param effect size estimates for PRT vs control on MBA scores -- hedges g
cats_delta = create_metadata_delta(cats);
delta_attr = cats_delta.mindbrain_attr_score_from_cats_delta;
clc

stats = mes(delta_attr(cats_delta.group==1), delta_attr(cats_delta.group==2), 'hedgesg')
stats = mes(delta_attr(cats_delta.group==1), delta_attr(cats_delta.group==3), 'hedgesg')
stats = mes(delta_attr(cats_delta.group==3), delta_attr(cats_delta.group==2), 'hedgesg')

%% Assocation of change in MBA scores and change in pain intensity  
clc
cats_delta = create_metadata_delta(cats);

cats_delta.pain_avg_delta_z = zscore(cats_delta.pain_avg_delta);
cats_delta.mindbrain_attr_score_from_cats_delta_z = zscore(cats_delta.mindbrain_attr_score_from_cats_delta);
cats_delta.grp_prt = cats_delta.group > 1;

fitglm(cats_delta, 'pain_avg_delta_z ~ mindbrain_attr_score_from_cats_delta_z * group + age + gender', 'CategoricalVars', 'group')

% follow-up: correlation of simple change scores: change in MBAS and pain within group or whole sample
clc, 
wh = cats_delta.group>0;
wh = cats_delta.group==1;

[r,p] = corr(cats_delta.pain_avg_delta(wh), delta_attr(wh))
[r,p] = corr(cats_delta.tsk11_delta(wh), delta_attr(wh))
[r,p] = corr(cats_delta.pcs_delta(wh), delta_attr(wh))

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
x = cats.mindbrain_attr_score_from_cats(cats.time==1);
mean(x), median(x), min(x), max(x)


%% assocation between gender and baseline attribution
% perm test since non-normal
wh = cats.time==1;
whg = cats.gender==1;

figure; 
histogram(cats.mindbrain_attr_score_from_cats(wh & whg))
hold on,
histogram(cats.mindbrain_attr_score_from_cats(wh & ~whg))

[p, observeddifference, effectsize] = permutationTest(cats.mindbrain_attr_score_from_cats(wh & whg), cats.mindbrain_attr_score_from_cats(wh & ~whg), 10000, 'plotresult', 0)

ttest2(cats.mindbrain_attr_score_from_cats(wh & whg), cats.mindbrain_attr_score_from_cats(wh & ~whg))

%% assocation between baseline attribution and X
wh = cats.time==1;
Y = cats.mindbrain_attr_score_from_cats(wh);


X = cats.backpain_length(wh);
%X = cats.age(wh);
X = cats.pain_avg(wh);
% 
% X = cats.tsk11(wh);
% X = cats.pcs(wh);
% 
% X = cats.sopa_emo(wh);

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

% just rename a variable as expected by the sub-function
cats_wide.mindbrain_attr_score_from_cats_baseline = cats_wide.mindbrain_attr_score_from_cats_t1_arm_1;

mediation_to_follow_up(1, cats_wide, 'mindbrain_attr_score_from_cats', 'pain_avg')


%% subfunctions

% input: array of categories
% output: 0/1 whether each category counts as mind/brain
function is_mb = cat_is_mindbrain(cat_array)
    
    is_mb = (cat_array == 'psychological') + (cat_array == 'brain') + (cat_array == 'stress');
    is_mb(isundefined(cat_array)) = NaN;
end    

function plot_pie(dat, fig_name)
    basedir = '/Users/yoni/Repositories/OLP4CBP/scripts/analyses/semantic_analyses';
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
    basedir = '/Users/yoni/Repositories/OLP4CBP/scripts/analyses/semantic_analyses';
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

% PRTvsTAU: 1 or 0 (0 does PRTvsOLP)
% mediation for M at post-tx, Outcome at all timepoints, controlling for
% baseline values of M and Outcome
function STATS =  mediation_to_follow_up(PRTvsTAU, outcomes_wide, Mname, Yname)

    outcome = get_outcome_data_from_final_acute_outcomes_wide_prt(outcomes_wide, Yname, 'noplot', 'nogpml', '1yearnoweekly', 'nosubtractbaseline');
    mymediator = get_outcome_data_from_final_acute_outcomes_wide_prt(outcomes_wide, Mname, 'noplot', 'nogpml', 'pre_post_only', 'nosubtractbaseline') ;  % substracts baseline by default
    posttx_ind = find(strcmp(outcome.xlab, ' Post-tx'));
    
    if PRTvsTAU     
        % PRT vs TAU
        X = [ones(size(outcome.datmat_prt, 1), 1); zeros(size(outcome.datmat_control, 1), 1)];
        outcome_datmat = [outcome.datmat_prt; outcome.datmat_control];
        mediator_datmat = [mymediator.datmat_prt; mymediator.datmat_control];
    else
        % PRT vs PLA
        X = [ones(size(outcome.datmat_prt, 1), 1); zeros(size(outcome.datmat_pla, 1), 1)];
        outcome_datmat = [outcome.datmat_prt; outcome.datmat_pla];
        mediator_datmat = [mymediator.datmat_prt; mymediator.datmat_pla];
    end
    
    M = mediator_datmat(:, posttx_ind); % mediator at post-tx
    cov = [outcome_datmat(:,1) mediator_datmat(:,1)]; % baseline vals of outcome and mediator
    
    clc, STATS = {};
    for i=1:6 % from post-tx to last timepiont

        print_header(['M = Delta ' Mname ' from baseline to post-tx'], ['Y = ' Yname ' at ' outcome.xlab{i + posttx_ind - 1}])
        Y = outcome_datmat(:, i + posttx_ind - 1); % outcome at given timepoint

        [~, stats] = mediation(X,scale(Y),scale(M), 'verbose',  'boot', 'bootsamples', 10000, 'covs', cov);
                 
        STATS.a(i, 1) = stats.z(1); % mean(1) ./ stats.ste(1);
        STATS.a_p(i, 1) = stats.p(1);

        STATS.b(i, 1) = stats.z(2); %./ stats.ste(2);
        STATS.b_p(i, 1) = stats.p(2);

        STATS.ab(i, 1) = stats.z(5); %./ stats.ste(5);
        STATS.ab_p(i, 1) = stats.p(5);

        STATS.perc_mediated(i, 1) = 100 .* (stats.mean(4) - stats.mean(3)) ./ stats.mean(4);
    end
end

% PRTvsTAU: 1 or 0 (0 does PRTvsOLP)
function STATS =  OLD_mediation_to_follow_up(PRTvsTAU, outcomes_wide, Mname, Yname, control_for_TSK11)

    outcome = get_outcome_data_from_final_acute_outcomes_wide_prt(outcomes_wide, Yname, 'noplot', 'nogpml', '1yearnoweekly', 'nosubtractbaseline');

    clc
    STATS = {};

    posttx_ind = find(strcmp(outcome.xlab, ' Post-tx'));
    
    mymediator = get_outcome_data_from_final_acute_outcomes_wide_prt(outcomes_wide, Mname, 'noplot', 'nogpml', 'pre_post_only', 'nosubtractbaseline') ;  % substracts baseline by default
 
    % include TSK11 as a coviarate
    cov = get_outcome_data_from_final_acute_outcomes_wide_prt(outcomes_wide, 'tsk11', 'noplot', 'nogpml', 'pre_post_only') ; % substracts baseline by default
    
    if PRTvsTAU     
        % PRT vs TAU
        X = [ones(size(outcome.datmat_prt, 1), 1); zeros(size(outcome.datmat_control, 1), 1)];
        outcome_datmat = [outcome.datmat_prt; outcome.datmat_control];
        mediator_datmat = [mymediator.datmat_prt; mymediator.datmat_control];
        cov_datmat = [cov.datmat_prt; cov.datmat_control];
    else
        % PRT vs PLA
        X = [ones(size(outcome.datmat_prt, 1), 1); zeros(size(outcome.datmat_pla, 1), 1)];
        outcome_datmat = [outcome.datmat_prt; outcome.datmat_pla];
        mediator_datmat = [mymediator.datmat_prt; mymediator.datmat_pla];
        cov_datmat = [cov.datmat_prt; cov.datmat_pla];
    end
    
    outcome_baseline = outcome_datmat(:,1); % used as covariate in all regressions below
    M = mediator_datmat(:, posttx_ind); % change from pre to post
    cov = cov_datmat(:, posttx_ind);
    
    for i=1:6 % from post-tx to last timepiont

        print_header(['M = Delta ' Mname ' from baseline to post-tx'], ['Y = ' Yname ' at ' outcome.xlab{i + posttx_ind - 1}])
        Y = outcome_datmat(:, i + posttx_ind - 1); % change from baseline to given timepoint

        if control_for_TSK11
            [~, stats] = mediation(X,scale(Y),scale(M), 'verbose',  'boot', 'bootsamples', 10000, 'covs', [outcome_baseline cov]);
        else
            [~, stats] = mediation(X,scale(Y),scale(M), 'verbose',  'boot', 'bootsamples', 10000, 'covs', outcome_baseline);
        end
            
        STATS.a(i, 1) = stats.z(1); % mean(1) ./ stats.ste(1);
        STATS.a_p(i, 1) = stats.p(1);

        STATS.b(i, 1) = stats.z(2); %./ stats.ste(2);
        STATS.b_p(i, 1) = stats.p(2);

        STATS.ab(i, 1) = stats.z(5); %./ stats.ste(5);
        STATS.ab_p(i, 1) = stats.p(5);

        STATS.perc_mediated(i, 1) = 100 .* (stats.mean(4) - stats.mean(3)) ./ stats.mean(4);
    end
end
