basedir = '/Users/yoni/Repositories/OLP4CBP/scripts/analyses/semantic_analyses';
figdir = fullfile(basedir, 'figures');
document_scores = readtable(fullfile(basedir, 'manual_codings', 'document_scores_from_parrot.csv'));
metadata = readtable(fullfile(basedir, 'manual_codings', 'IPQlastitem_final MERGED LONG FORMAT T2 only'));

%% "merge" with metadata. confirmed that they are lined up
head(document_scores)
head(metadata)

document_scores.id = metadata.id;
document_scores.group = metadata.group;

%% plot

placebo_color = [.5 0 .5];
tau_color = [.3 .3 .3];
prt_color = [0 0 1];

colormat = repmat(tau_color, height(document_scores), 1); % set all to black
colormat(document_scores.group==2, :) = repmat(placebo_color, sum(document_scores.group==2), 1);
colormat(document_scores.group==1, :) = repmat(prt_color, sum(document_scores.group==1), 1); 

%%
create_figure('myfig')
scatter(document_scores.X0, document_scores.X1, 88, colormat, 'filled')



