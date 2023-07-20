clc
basedir = '/Users/yoni/Repositories/OLP4CBP/scripts/analyses/semantic_analyses/manual_codings';
attr = readtable(fullfile(basedir, 'IPQlastitem_final FINAL CODED.csv'));
attr.rowkey = [];
attr.Properties.VariableNames{1} = 'id';
attr.Properties.VariableNames(3:8) = {'attributed_cause_1', 'attributed_cause_2', 'attributed_cause_3', 'attributed_code_1', 'attributed_code_2', 'attributed_code_3'};
attr.redcap_event_name = categorical(attr.redcap_event_name);

head(attr)

% merge CATEGORIES with metadata
cats = readtable(fullfile(basedir, 'IPQlastitem categories.csv'));
cats.rowkey = [];
cats.Properties.VariableNames{1} = 'id';
cats.Properties.VariableNames(3:8) = {'attributed_cause_1', 'attributed_cause_2', 'attributed_cause_3', 'cat1', 'cat2', 'cat3'};
cats.redcap_event_name = categorical(cats.redcap_event_name);


% load in outcomes data
basedir = '/Users/yoni/Repositories/OLP4CBP/data';
load(fullfile(basedir, 'final_12mo_outcomes_long.mat'))
load(fullfile(basedir, 'final_12mo_outcomes_wide.mat'))
load(fullfile(basedir, 'all_subjects_outcomes_demographics.mat'))

demographics = all_subjects_outcomes_demographics(all_subjects_outcomes_demographics.time==1 & all_subjects_outcomes_demographics.is_patient,:);
demographics = demographics(:, [1 22:end]);

% a number of Ps completed the IPQ but then were not randomized (eg, fMRI
% problem, claustrophobia, etc.). Drop these non randomized Ps
idstodrop = setdiff(attr.id, outcomes_wide.id);
attr(ismember(attr.id, idstodrop),:) = []; 
cats(ismember(cats.id, idstodrop),:) = []; 


%% merge with outcomes_long

% this merge drops all timepoints in outcomes_long besides T1 and T2
attr_long = join(attr, outcomes_long,  'Keys', {'id' 'redcap_event_name'});
attr_long = join(attr_long, demographics );
head(attr_long)

%% merge cat with outcomes_long
cats_long = join(cats, outcomes_long,  'Keys', {'id' 'redcap_event_name'});
cats_long = join(cats_long, demographics );
head(cats_long)

%% merge with outcomes_wide 
attr_wide = unstack(attr, attr.Properties.VariableNames(3:end), 'redcap_event_name');
attr_wide = join(attr_wide, outcomes_wide, 'Keys', {'id'});
attr_wide = join(attr_wide, demographics);

cats_wide = unstack(cats, cats.Properties.VariableNames(3:end), 'redcap_event_name');
cats_wide = join(cats_wide, outcomes_wide, 'Keys', {'id'});
cats_wide = join(cats_wide, demographics);

%% save data files

basedir = '/Users/yoni/Repositories/OLP4CBP/scripts/analyses/semantic_analyses/manual_codings';
writetable(attr_wide, fullfile(basedir, 'IPQlastitem_final MERGED WIDE FORMAT.csv'));
writetable(attr_long, fullfile(basedir, 'IPQlastitem_final MERGED LONG FORMAT.csv'));

writetable(cats_long, fullfile(basedir, 'IPQlastitem_final MERGED LONG CATEGORIES.csv'));
writetable(cats_wide, fullfile(basedir, 'IPQlastitem_final MERGED WIDE CATEGORIES.csv'));

