basedir = '/Users/yoni/Repositories/OLP4CBP/scripts/analyses/semantic_analyses';
figdir = fullfile(basedir, 'figures');
attr_wide = readtable(fullfile(basedir, 'manual_codings', 'IPQlastitem_final MERGED WIDE FORMAT.csv'));
attr_long = readtable(fullfile(basedir, 'manual_codings', 'IPQlastitem_final MERGED LONG FORMAT.csv'));

% code for time, in long form
a=categorical(attr_long.redcap_event_name);
attr_long.time = grp2idx(a);

%% create an MB score across the top 3 reasons
attr_wide.attr_mb_score_t1 = strcmp(attr_wide.attributed_code_1_t1_arm_1,'MB') + strcmp(attr_wide.attributed_code_2_t1_arm_1,'MB') + strcmp(attr_wide.attributed_code_3_t1_arm_1,'MB');
attr_wide.attr_mb_score_t2 = strcmp(attr_wide.attributed_code_1_t2_arm_1,'MB') + strcmp(attr_wide.attributed_code_2_t2_arm_1,'MB') + strcmp(attr_wide.attributed_code_3_t2_arm_1,'MB');

attr_long.attr_mb_score = strcmp(attr_long.attributed_code_1,'MB') + strcmp(attr_long.attributed_code_2,'MB') + strcmp(attr_long.attributed_code_3,'MB');

%% distribution of MB scores at T1 by reason number. 
% Slightly LESS MB ratings for code 1 relative to codes 2 and 3 at T1
% At T2, pretty even distribution across codes

A = categorical(attr_wide.attributed_code_1_t2_arm_1, {'<undefined>' 'STR' 'MB'}, {'N/A' 'STR' 'MB'} );
clc
unique(A)
countcats(A)

create_figure('basics', 1,2)
histogram(attr_wide.attr_mb_score_t1);
hold on,
histogram(attr_wide.attr_mb_score_t2);
legend({'T1' 'T2'})

subplot(1,2,2)
histogram(attr_wide.attr_mb_score_t2 - attr_wide.attr_mb_score_t1)
title('shift per person')

%% assocation between baseline attribution and X
X = attr_wide.backpain_length;
%X = attr_wide.age;
%X = attr_wide.pain_avg_baseline;

%X = attr_wide.tsk11_baseline;
%X = attr_wide.pcs_baseline;

%X = attr_wide.sopa_emo_t1_arm_1;

clc
if 0
    figure, scatterhist(X, attr_wide.attr_mb_score_t1);
    set(gca, 'FontSize', 24), ylabel('Baseline attribution');xlabel('X'); lsline
end
[B,STATS] = robustfit(zscore(attr_wide.attr_mb_score_t1), zscore(X));
B(2), STATS.p(2)

[r,p]=corr(X, attr_wide.attr_mb_score_t1, 'Type', 'Spearman')

%% assocation between post-tx attribution and X in PRT group

wh = attr_wide.group==1;
X = attr_wide.tsk11_t2_arm_1;
%X = attr_wide.sopa_emo_t2_arm_1;

clc
figure, scatterhist(X(wh), attr_wide.attr_mb_score_t2(wh));
set(gca, 'FontSize', 24), ylabel('Baseline attribution');xlabel('X'); lsline
[B,STATS] = robustfit(scale(attr_wide.attr_mb_score_t2(wh)), scale(X(wh)));
B(2), STATS.p(2)

%% assocation between gender and baseline attribution
% perm test since non-normal
wh = attr_wide.gender==1;
figure; violinplot({attr_wide.attr_mb_score_t1(wh), attr_wide.attr_mb_score_t1(~wh)} )
[p, observeddifference, effectsize] = permutationTest(attr_wide.attr_mb_score_t1(wh), attr_wide.attr_mb_score_t1(~wh), 10000, 'plotresult', 1)


%% non-param effect size estimates for PRT vs control

delta_attr = attr_wide.attr_mb_score_t2 - attr_wide.attr_mb_score_t1;
clc
disp('PRT vs PLA')
[p, observeddifference, effectsize] = permutationTest(delta_attr(attr_wide.group==1), delta_attr(attr_wide.group==2), 10000, 'plotresult', 1)
disp('PRT vs UC')
[p, observeddifference, effectsize] = permutationTest(delta_attr(attr_wide.group==1), delta_attr(attr_wide.group==3), 10000)

disp('PLA vs UC')
[p, observeddifference, effectsize] = permutationTest(delta_attr(attr_wide.group==2), delta_attr(attr_wide.group==3), 10000)


%% pre-to-post-tx changes in attributions by X in PRT group

wh = attr_wide.group==1;
delta_tsk11 = attr_wide.tsk11_t2_arm_1 - attr_wide.tsk11_baseline;
delta_pain = attr_wide.pain_avg_t2_arm_1 - attr_wide.pain_avg_baseline;
delta_odi = attr_wide.odi_t2_arm_1 - attr_wide.odi_baseline;

X = delta_tsk11;

if 0
    figure, scatterhist(X, delta_attr);
    set(gca, 'FontSize', 24), ylabel('Delta attribution');xlabel('Delta X'); lsline
end

clc
[B,STATS] = robustfit(scale(delta_attr(wh)), scale(X(wh)));
B(2), STATS.p(2)

[r,p] = corr(delta_attr(wh), X(wh), 'Type', 'Spearman', 'rows', 'complete')


%% multiple regression: delta pain = delta MBAscore + delta TSK11
clc
[B,STATS] = robustfit(scale([delta_attr(wh) delta_tsk11(wh)]), scale(delta_pain(wh)));
B
STATS.p

mdl = fitglm(scale([delta_attr(wh) delta_tsk11(wh)]), scale(delta_pain(wh)))


%% mediation models
attr_wide.attr_mb_score_baseline = attr_wide.attr_mb_score_t1;
attr_wide.attr_mb_score_t2_arm_1 = attr_wide.attr_mb_score_t2;

mediation_to_follow_up(1, attr_wide, 'attr_mb_score', 'pain_avg', 1)

% sanity check:
% mediation_to_follow_up(1, attr_wide, 'tsk11', 'pain_avg')

%% test whether 1yr follow-up results hold when control for TSK11
% TSK11 is a stronger predictor -- but is it sig stronger? we don't know
fitglm([attr_wide.pain_avg_baseline delta_attr delta_tsk11], attr_wide.pain_avg_x12_month_follow_up_arm_1, 'VarNames', {'baseline pain', 'delta MBA', 'delta TSK', 'pain 1yr'}, 'Exclude', attr_wide.group~=1)


%% ML model 

model_olp4cbp_outcomes(attr_long.attr_mb_score, attr_long, 'standardize')

% younger ppl have more mind body p = .06
% female have higher mind body,  p = .06  (M=1, F=2)

%% Plot the group by time
[h1,h2,h3] = plot_olp4cbp_group_by_time(attr_long, 'attr_mb_score');

figure(h1)
set(h1, 'Position', [440   586   290   225])

ylabel('Mind-brain attribution (0-3)')
legend('Location', 'NW')
print_pdf(h1, fullfile(figdir, 'attributions_grp_by_time.pdf'))

close(h3)

figure(h2)
set(gca, 'FontSize', 16)
xlabel('Δ Pain, Pre-to-post-tx');
ylabel({'Mind-brain attribution', 'Δ Pre-to-post-tx'});
%legend({'PRT' 'Placebo' 'No Tx'});
legend off
set(gca, 'FontSize', 24);
print_pdf(h2, fullfile(figdir, 'attributions_ind_diffs.pdf'))

%% print words to file for Word Cloud

f1 = fopen(fullfile(figdir, 'PRT_T1_words.txt'), 'w');
sa = attr_long.attributed_cause_1(attr_long.group==1 & attr_long.time==1);
sb = attr_long.attributed_cause_2(attr_long.group==1 & attr_long.time==1);
sc = attr_long.attributed_cause_3(attr_long.group==1 & attr_long.time==1);
s = cat(1, sa, sb, sc);
fprintf(f1, '%s\n', s{:});

f1 = fopen(fullfile(figdir, 'PRT_T2_words.txt'), 'w');
sa = attr_long.attributed_cause_1(attr_long.group==1 & attr_long.time==2);
sb = attr_long.attributed_cause_2(attr_long.group==1 & attr_long.time==2);
sc = attr_long.attributed_cause_3(attr_long.group==1 & attr_long.time==2);
s = cat(1, sa, sb, sc);
fprintf(f1, '%s\n', s{:});

f1 = fopen(fullfile(figdir, 'Control_T1_words.txt'), 'w');
sa = attr_long.attributed_cause_1(attr_long.group~=1 & attr_long.time==1);
sb = attr_long.attributed_cause_2(attr_long.group~=1 & attr_long.time==1);
sc = attr_long.attributed_cause_3(attr_long.group~=1 & attr_long.time==1);
s = cat(1, sa, sb, sc);
fprintf(f1, '%s\n', s{:});

f1 = fopen(fullfile(figdir, 'Control_T2_words.txt'), 'w');
sa = attr_long.attributed_cause_1(attr_long.group~=1 & attr_long.time==2);
sb = attr_long.attributed_cause_2(attr_long.group~=1 & attr_long.time==2);
sc = attr_long.attributed_cause_3(attr_long.group~=1 & attr_long.time==2);
s = cat(1, sa, sb, sc);
fprintf(f1, '%s\n', s{:});

%% SUBFUNCTIONS

% PRTvsTAU: 1 or 0 (0 does PRTvsOLP)
function STATS =  mediation_to_follow_up(PRTvsTAU, outcomes_wide, Mname, Yname, control_for_TSK11)

    outcome = get_outcome_data_from_final_acute_outcomes_wide_prt(outcomes_wide, Yname, 'noplot', 'nogpml', '1yearnoweekly', 'nosubtractbaseline');

    clc
    STATS = {};

    posttx_ind = find(strcmp(outcome.xlab, ' Post-tx'));
    
    mymediator = get_outcome_data_from_final_acute_outcomes_wide_prt(outcomes_wide, Mname, 'noplot', 'nogpml', 'pre_post_only') ;  % substracts baseline by default
 
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
