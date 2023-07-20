% This script counts how many expert list words are present in attributions
% Should do a dictionary expansion to include SYNONYMS and RELATED WORDS
% This works in a training subset, balanced on group and time

%% load attributions

basedir = '/Users/yoni/Repositories/OLP4CBP/scripts/analyses/semantic_analyses';
figdir = fullfile(basedir, 'figures');
attr_long = readtable(fullfile(basedir, 'manual_codings', 'IPQlastitem_final MERGED LONG FORMAT.csv'));

% code for time, in long form
a=categorical(attr_long.redcap_event_name);
attr_long.time = grp2idx(a);

% load preproc'd expert list
expert_list = readtable(fullfile(basedir, 'expert_lists', 'mindbrain_list_preprocessed_cleaned.csv'), 'ReadVariableNames', false);
expert_list = expert_list.Var1;


%% word length of documents
tmp = [tokenizedDocument(attr_long.attributed_cause_1); tokenizedDocument(attr_long.attributed_cause_2); tokenizedDocument(attr_long.attributed_cause_3)];
tmp2 = doclength(tmp);
clc
mean(tmp2), min(tmp2), max(tmp2), std(tmp2), median(tmp2)

%% preproc attributions and then compare to expert list

attr_long.attr1 = joinWords(preprocessTextData(attr_long.attributed_cause_1, 1));
attr_long.attr2 = joinWords(preprocessTextData(attr_long.attributed_cause_2, 1));
attr_long.attr3 = joinWords(preprocessTextData(attr_long.attributed_cause_3, 1));

% manual check for correct token counting -- #OK
% ism = contains(attr_long.attr1, expert_list);
% tmp = table(ism, attr_long.attr1);
% tmp(133:150,:)

%% compute mind brain attribution scores

% warning: contains('diet','die') == 1

% token count at each observation -- 0 to 3
attr_long.mb_token_count = contains(attr_long.attr1, expert_list) + contains(attr_long.attr2, expert_list) + contains(attr_long.attr3, expert_list);


%% compare token-counting to manual scores from categories -- Cohen's kappa

manual_scores = readtable(fullfile(basedir, 'manual_codings', 'IPQlastitem_final MERGED LONG CATEGORIES.csv'));

% confirmed that they do line up

a = manual_scores.mindbrain_attr_score_from_cats;
b = attr_long.mb_token_count;
time1 = attr_long.time==1;

% Cohen's kappa, from: http://www.mathworks.com/matlabcentral/fileexchange/15365
[cm, order] = confusionmat(a(~time1),b(~time1));
kappa(cm, 1) % use linear weighting, as categories are ordinal not nominal

figure
confusionchart(cm, order)


%% which expert list words have the most hits (see also analyses_of_words.m)


% find the topkwords in the attributions and the freq count of each
% find the ones that are in the expert list

% PRT baseline
wh = attr_long.group == 1 & strcmp(attr_long.redcap_event_name, 't1_arm_1');
mywords = [attr_long.attributed_cause_1(wh); attr_long.attributed_cause_2(wh); attr_long.attributed_cause_3(wh)];
mytopwordsT1 = get_count_for_top_expert_derived_words(mywords, expert_list)
mytopwordsT1.Properties.VariableNames{2} = 'Count_T1';

% PRT post-tx
wh = attr_long.group == 1 & strcmp(attr_long.redcap_event_name, 't2_arm_1');
mywords = [attr_long.attributed_cause_1(wh); attr_long.attributed_cause_2(wh); attr_long.attributed_cause_3(wh)];
mytopwordsT2 = get_count_for_top_expert_derived_words(mywords, expert_list)
mytopwordsT2.Properties.VariableNames{2} = 'Count_T2';

% join T1 and T2
joined = outerjoin(mytopwordsT1, mytopwordsT2);
joined.Count_T1(isnan(joined.Count_T1)) = 0;
joined.Count_T2(isnan(joined.Count_T2)) = 0;
joined.Word_mytopwordsT1(ismissing(joined.Word_mytopwordsT1)) = joined.Word_mytopwordsT2(ismissing(joined.Word_mytopwordsT1));
joined.Word_mytopwordsT2 = [];
joined

%% plot T1 and T2 stacked/together -- doesn't look good

create_figure('top expert words'), bar(joined{:, 2:3}, 'LineWidth', 1)
set(gca, 'XTick', 1:height(joined), 'XTickLabel', joined.Word_mytopwordsT1, 'FontSize', 16)
ylabel('# of attributions'), ylim([0 20])
legend({'PRT Baseline', 'PRT Post-tx'}, 'Location', 'best')

% update fname when
fig_fname = fullfile(figdir, 'top_expert_words_hits_T1_T2.pdf');
saveas(gcf, fig_fname);
%print(gcf, '-dpdf', fig_fname);


%% plot T1 and T2 sequentially with a gap in the middle

mytopwordsT1.Properties.VariableNames{2} = 'Count';
mytopwordsT2.Properties.VariableNames{2} = 'Count';
blank = table();
blank.Word(1) = ' '; blank.Word(2)=  ' ';
blank.Count = [0 0]';

sequential = [mytopwordsT1; blank; mytopwordsT2];
%%
create_figure('top expert words T1'), 
h = bar(1:5, mytopwordsT1.Count, 'FaceColor', [.6 .6 .8]);
h = bar(8:17, mytopwordsT2.Count, 'FaceColor', [.4 .4 .7]);

mylabels = strcat('''', sequential.Word, '''');
mylabels(6) = ''; mylabels(7) = '';

set(gca, 'XTick', 1:height(sequential), 'XTickLabel', mylabels, 'FontSize', 18, 'TickLength',[0 0])
ylabel('# of attributions'), ylim([0 18])
legend({'PRT Pre-tx', 'PRT Post-tx'}, 'Location', 'best')
set(gcf, 'Position', [       1105         135         633         275]);

% update fname when
fig_fname = fullfile(figdir, 'top_expert_words_hits_T1_T2_seq.pdf');
%saveas(gcf, fig_fname);
print(gcf, '-dpdf', fig_fname);



%% repeat the cell above for attr_wide, so have it in that format too

attr_wide = readtable(fullfile(basedir, 'manual_codings', 'IPQlastitem_final MERGED WIDE FORMAT.csv'));
fields = {'attributed_cause_1_t1_arm_1' 'attributed_cause_2_t1_arm_1' 'attributed_cause_3_t1_arm_1' 'attributed_cause_1_t2_arm_1' 'attributed_cause_2_t2_arm_1' 'attributed_cause_3_t2_arm_1'};

% preprocess
for f = fields
    attr_wide.([f{1} '_raw']) = attr_wide.(f{1});
    attr_wide.(f{1}) = joinWords(preprocessTextData(attr_wide.(f{1}), 1));
    
    % score whether this one is MB or not
    attr_wide.([f{1} '_is_mb']) = contains(attr_wide.(f{1}), expert_list);
    
end

% score each timepoint. Add in NaN for missing data
attr_wide.mb_token_count_baseline = contains(attr_wide.(fields{1}), expert_list) + contains(attr_wide.(fields{2}), expert_list) + contains(attr_wide.(fields{3}), expert_list);
attr_wide.mb_token_count_t2_arm_1 = contains(attr_wide.(fields{4}), expert_list) + contains(attr_wide.(fields{5}), expert_list) + contains(attr_wide.(fields{6}), expert_list);
attr_wide.mb_token_count_t2_arm_1(isnan(attr_wide.pain_avg_t2_arm_1)) = NaN;



%% baseline distribution of MB scores

histc(attr_wide.mb_token_count_baseline, 0:3)
histc(attr_wide.mb_token_count_baseline, 0:3) ./ height(attr_wide) * 100

%% correlation of NLP and manual scores

% create an MB score across the top 3 reasons from the manual codings
attr_wide.attr_mb_score_t1 = strcmp(attr_wide.attributed_code_1_t1_arm_1,'MB') + strcmp(attr_wide.attributed_code_2_t1_arm_1,'MB') + strcmp(attr_wide.attributed_code_3_t1_arm_1,'MB');
attr_wide.attr_mb_score_t2 = strcmp(attr_wide.attributed_code_1_t2_arm_1,'MB') + strcmp(attr_wide.attributed_code_2_t2_arm_1,'MB') + strcmp(attr_wide.attributed_code_3_t2_arm_1,'MB');

% correlate with NLP codings
[r,p] = corr(attr_wide.mb_token_count_baseline, attr_wide.attr_mb_score_t1, 'Type', 'Spearman')
[r,p] = corr(attr_wide.mb_token_count_t2_arm_1, attr_wide.attr_mb_score_t2, 'Type', 'Spearman', 'rows', 'complete')

%% assocation between gender and baseline attribution
% perm test since non-normal
wh = attr_wide.gender==1;
figure; violinplot({attr_wide.mb_token_count_baseline(wh), attr_wide.mb_token_count_baseline(~wh)} )
[p, observeddifference, effectsize] = permutationTest(attr_wide.mb_token_count_baseline(wh), attr_wide.mb_token_count_baseline(~wh), 10000, 'plotresult', 1)

%% assocation between baseline attribution and X
X = attr_wide.backpain_length;
X = attr_wide.age;
%X = attr_wide.pain_avg_baseline;

% X = attr_wide.tsk11_baseline;
% X = attr_wide.pcs_baseline;
% X = attr_wide.sopa_emo_t1_arm_1;

clc
if 0
    figure, scatterhist(X, attr_wide.mb_token_count_baseline);
    set(gca, 'FontSize', 24), ylabel('Baseline attribution');xlabel('X'); lsline
end
[B,STATS] = robustfit(zscore(attr_wide.mb_token_count_baseline), zscore(X));
B(2), STATS.p(2)

%[r,p]=corr(X, attr_wide.mb_token_count_baseline, 'Type', 'Spearman')

%% model
attr_long.grp = attr_long.group;
model_olp4cbp_outcomes(attr_long.mb_token_count, attr_long)

%% plot
[h1,h2,h3, tabl_delta] = plot_olp4cbp_group_by_time(attr_long, 'mb_token_count');

figure(h1)
set(h1, 'Position', [440   586   290   225])

ylabel('Mind-brain attribution (0-3)')
legend('Location', 'NW')
print_pdf(h1, fullfile(figdir, 'attributions_expert-derived_grp_by_time.pdf'))

figure(h2)
set(gca, 'FontSize', 16)
xlabel('Δ Pain, Pre-to-post-tx');
ylabel({'Mind-brain attribution', 'Δ Pre-to-post-tx'});
%legend({'PRT' 'Placebo' 'No Tx'});
legend off
set(gca, 'FontSize', 24);
print_pdf(h2, fullfile(figdir, 'attributions_ind_diffs.pdf'))

%% effect sizes

% drop missing
attr_wide_completers = attr_wide(~isnan(attr_wide.mb_token_count_t2_arm_1),:);

delta_attr = attr_wide_completers.mb_token_count_t2_arm_1 - attr_wide_completers.mb_token_count_baseline;
clc
disp('PRT vs PLA')
[p, observeddifference, effectsize] = permutationTest(delta_attr(attr_wide_completers.group==1), delta_attr(attr_wide_completers.group==2), 10000)
disp('PRT vs UC')
[p, observeddifference, effectsize] = permutationTest(delta_attr(attr_wide_completers.group==1), delta_attr(attr_wide_completers.group==3), 10000)

disp('PLA vs UC')
[p, observeddifference, effectsize] = permutationTest(delta_attr(attr_wide_completers.group==2), delta_attr(attr_wide_completers.group==3), 10000)


%% pre-to-post-tx changes in attributions by X in PRT group

wh = attr_wide_completers.group==1;
delta_tsk11 = attr_wide_completers.tsk11_t2_arm_1 - attr_wide_completers.tsk11_baseline;
delta_pain = attr_wide_completers.pain_avg_t2_arm_1 - attr_wide_completers.pain_avg_baseline;
delta_odi = attr_wide_completers.odi_t2_arm_1 - attr_wide_completers.odi_baseline;

X = delta_pain;

clc
if 0
    figure, scatterhist(X, delta_attr);
    set(gca, 'FontSize', 24), ylabel('Delta attribution');xlabel('Delta X'); lsline
end

% robust fit delta X delta Y
[B,STATS] = robustfit(scale(delta_attr(wh)), scale(X(wh)));
B(2), STATS.p(2)

% corr delta X delta Y
[r,p] = corr(delta_attr(wh), X(wh), 'Type', 'Spearman', 'rows', 'complete')

%% Time 2 pain ~ Time 1 pain + T1 attr + T2 attr

%[B,STATS] = robustfit(attr_wide_completers{:, {'pain_avg_t1_arm_1' 'mb_token_count_baseline' 'mb_token_count_t2_arm_1'}}, attr_wide_completers.pain_avg_t2_arm_1)
clc
attr_wide_completers.delta_attr = attr_wide_completers.mb_token_count_t2_arm_1 - attr_wide_completers.mb_token_count_baseline;
attr_wide_completers.PRT_vs_control = attr_wide_completers.group==1;
attr_wide_completers.pla_vs_uc = double(attr_wide_completers.group==2) + (attr_wide_completers.group==3)*-1;

% for each additional attribution scored as MB at T2, 0.75 (of 10) units of less
% pain at T2
fitglm(attr_wide_completers, 'pain_avg_t2_arm_1 ~ pain_avg_t1_arm_1 + mb_token_count_baseline + mb_token_count_t2_arm_1', 'Exclude', ~wh)
fitglm(attr_wide_completers, 'pain_avg_t2_arm_1 ~ pain_avg_t1_arm_1 + mb_token_count_baseline + mb_token_count_t2_arm_1 + PRT_vs_control ')


attr_wide_completers.delta_attr = attr_wide_completers.mb_token_count_t2_arm_1 - attr_wide_completers.mb_token_count_baseline;
fitglm(attr_wide_completers, 'pain_avg_t2_arm_1 ~ pain_avg_t1_arm_1 + delta_attr')

attr_wide_completers.delta_pain = attr_wide_completers.pain_avg_t2_arm_1 - attr_wide_completers.pain_avg_baseline;
fitglm(attr_wide_completers, 'delta_pain ~  delta_attr + PRT_vs_control + PRT_vs_control*delta_attr')


%% multiple regression: delta pain = delta MBAscore + delta TSK11
clc
[B,STATS] = robustfit(scale([delta_attr(wh) delta_tsk11(wh)]), scale(delta_pain(wh)));
B
STATS.p

mdl = fitglm(scale([delta_attr(wh) delta_tsk11(wh)]), scale(delta_pain(wh)))


%% PF-NPF -- no diffs when appropriately drop Ss who withdrew. It might work with manual scoring tho.
wh = attr_wide_completers.group==1;
pfnpf = attr_wide_completers.pain_avg_t2_arm_1 <= 1;
pf = attr_wide_completers.pain_avg_t2_arm_1 < 1;
inpain = attr_wide_completers.pain_avg_t2_arm_1 >= 3;

Y = attr_wide_completers.mb_token_count_baseline;

figure;
violinplot({Y(wh & pf), Y(wh & ~pf)});

[p, ~, stats] = ranksum(Y(wh & pf), Y(wh & ~pf))

attr_wide_completers.pf = pf;
fitglm(attr_wide_completers, 'mb_token_count_t2_arm_1 ~ mb_token_count_baseline + pf', 'Exclude', ~wh)

%% look at the raw words for PF at T1 and T2

painfree_tabl = table();
painfree_tabl.t1 = attr_wide_completers.attributed_cause_1_t1_arm_1_raw(wh & pf);
painfree_tabl.t2 = attr_wide_completers.attributed_cause_1_t2_arm_1_raw(wh & pf);
painfree_tabl

%% look at the raw words for PF vs in pain at T2
clc
disp('PRT group, pain >=3, attribution 3 at post-tx:')
disp(attr_wide_completers.attributed_cause_3_t2_arm_1_raw(wh & inpain))

disp('PRT group, pain = 0, attribution 3 at post-tx:')
disp(attr_wide_completers.attributed_cause_3_t2_arm_1_raw(wh & pf))


%% look at attributions 1 - 3 separately

% NLP
sum(attr_wide_completers{wh & pf, {'attributed_cause_1_t1_arm_1_is_mb' 'attributed_cause_2_t1_arm_1_is_mb' 'attributed_cause_3_t1_arm_1_is_mb'}})
sum(attr_wide_completers{wh & pf, {'attributed_cause_1_t2_arm_1_is_mb' 'attributed_cause_2_t2_arm_1_is_mb' 'attributed_cause_3_t2_arm_1_is_mb'}})

% manual
a = attr_wide_completers(wh & pf,:);
sum([ strcmp(a.attributed_code_1_t1_arm_1,'MB') strcmp(a.attributed_code_2_t1_arm_1,'MB') strcmp(a.attributed_code_3_t1_arm_1,'MB')])
sum([ strcmp(a.attributed_code_1_t2_arm_1,'MB') strcmp(a.attributed_code_2_t2_arm_1,'MB') strcmp(a.attributed_code_3_t2_arm_1,'MB')])


%%
figure
beeswarm(pfnpf(wh), attr_wide.mb_token_count_t2_arm_1(wh))

%% test mediation

mediation_to_follow_up(1, attr_wide, 'mb_token_count', 'pain_avg', 0)


%% SUBFUNCTIONS

% input: cell array of a bunch of words
% output: table with the top words from the expert derived list found in
% the input words, plus frequency count. Thresholded to return words with
% frequency >= 2
%
% NOTE: ignores any structure to the data, eg, ignores participants, which
% attribution, etc.
function tabl = get_count_for_top_expert_derived_words(mywords, expert_list)

    mybag = bagOfWords(preprocessTextData(mywords, 1));
    topwords_table = topkwords(mybag, 80);

    % awk hack to find top words that are in the expert list. can't use
    % contains because   contains('diet','die') == 1 . so I am using strcmp,
    % but then need to loop
    is_in_expert_list = false();
    for i=1:height(topwords_table)
        is_in_expert_list(i) = logical(sum(strcmp(topwords_table.Word(i), expert_list)));
    end

    tabl = topwords_table(is_in_expert_list,:);
    
    tabl = tabl(tabl.Count > 1,:);

end

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
