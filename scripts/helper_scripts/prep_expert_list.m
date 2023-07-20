basedir = '/Users/yoni/Repositories/OLP4CBP/scripts/analyses/semantic_analyses';
tabl = readtable(fullfile(basedir, 'expert_lists', 'Mindbrain_words.csv'), 'ReadVariableNames', false);

% remove duplicate words. I manually removed duplicate lemmas that remain
% after dup words removed
[a,b]=unique(lower(tabl.Var1));
nodupwords = tabl.Var1(b);

% lemmatize etc
clc
list = preprocessTextData(nodupwords)

writeTextDocument(list, fullfile(basedir, 'expert_lists', 'mindbrain_list_preprocessed.csv'))

% after this: undo wrong spelling corrections (eg maintain OCD, PPD), and
% remove a few that are not mind brain attributions (eg rain, weather) and
% save as 'mindbrain_list_preprocessed_cleaned.csv'