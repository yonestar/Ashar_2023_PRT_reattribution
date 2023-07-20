basedir = '/Users/yoni/Repositories/OLP4CBP/scripts/analyses/semantic_analyses/';

tabl = readtable(fullfile(basedir, 'parrot', 'word_scores_1st_dim_from_parrot.csv'));

%% thresh on most common words

mytabl = tabl(tabl.scores_word_counts > 4, :);

%% plot line

create_figure('word line')
scatter(mytabl.scores_word_scores___2_, zeros(1, height(mytabl)))
labelpoints(mytabl.scores_word_scores___2_, zeros(1, height(mytabl)), mytabl.scores_vocab, 'rotation', 50, 'FontSize', 24)

%% plot in ranked order

sorted = sortrows(mytabl, 2);
sorted.rank = (1:height(sorted))'

create_figure('word line')
scatter(sorted.rank, zeros(1, height(sorted)), '.w')
labelpoints(sorted.rank, zeros(1, height(sorted)), sorted.scores_vocab, 'rotation', 50, 'FontSize', 16, 'buffer', .7, 'position', 'N')
axis off

set(gcf, 'Position',[ 225       1147         816          83])

%% save fig

saveas(gcf, fullfile(basedir, 'figures', 'word_line.pdf'))