repodir = '/Users/yoni/Repositories/OLP4CBP/';
basedir = fullfile(repodir, 'scripts', 'analyses', 'semantic_analyses');

tabl = readtable(fullfile(basedir, 'IPQlastitem_final.csv'));
head(tabl)

%% add a unique blind identifier to each row, for later re-merging

tabl.rowkey = (1:height(tabl))';
head(tabl)

% save it back to file with the key
writetable(tabl, fullfile(basedir, 'IPQlastitem_final.csv'));

%% remove the identifying info, permute, and save it out

tabl_deid = tabl(:,[6 3:5]);
tabl_deid = tabl_deid( randperm(height(tabl)), :);
head(tabl_deid)
%%
writetable(tabl_deid, fullfile(basedir, 'IPQlastitem_final_deid.csv'));
