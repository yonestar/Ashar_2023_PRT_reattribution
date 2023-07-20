basedir = '/Users/yoni/Repositories/OLP4CBP/scripts/analyses/semantic_analyses/manual_codings';
key = readtable(fullfile(basedir, 'IPQlastitem_final.csv'));
head(key)

deid = readtable(fullfile(basedir, 'IPQlastitem_final_deid FINAL.csv'));

reid = join(key, deid);

head(reid)

%% save out re-identified file

writetable(reid, fullfile(basedir, 'IPQlastitem_final FINAL CODED.csv'))



%% repeat for the categories data

deid = readtable(fullfile(basedir, 'IPQlastitem_final_deid categories reconciled.csv'));
reid = join(key, deid);
clc
head(deid)
head(reid)

writetable(reid, fullfile(basedir, 'IPQlastitem categories.csv'))
