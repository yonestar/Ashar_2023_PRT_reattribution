clear all

basedir = '/Users/yoni/Repositories/Ashar_2023_PRT_reattribution/';

tabl = readtable(fullfile(basedir, 'scripts', 'parrot', 'document_scores_from_parrot.csv'));

% import this for metadata used below
cats_wide = readtable(fullfile(basedir, 'data', 'manual_codings', 'IPQlastitem_final MERGED WIDE CATEGORIES.csv'));
cats_wide = cats_wide(~isnan(cats_wide.pain_avg_t2_arm_1),:); % only keep completers

mymeta = cats_wide(:,{'id' 'group' 'pain_avg_baseline' 'pain_avg_t2_arm_1' 'age' 'gender'});
tabl.id = cats_wide.id;

% merge in metadata with text scaling data
tabl = join(tabl, mymeta);


%% scatter in first dim

create_figure('myscat')

jitterAmount = 0.5;
myzeros = zeros(1, height(tabl));
jitterValuesY = 2*(rand(size(myzeros))-0.5)*jitterAmount;   % +/-jitterAmount max

placebo_color = [.5 0 .5];
tau_color = [.3 .3 .3];
colors = [0 0 .8; placebo_color; tau_color] + repmat([.2 .2 .2], 3, 1);

h = gscatter(tabl.X1, jitterValuesY, tabl.group, zeros(3, 3), 'o', 14);
for n = 1:length(h)
  set(h(n), 'MarkerFaceColor', colors(n,:), 'LineWidth', 1.2);
end
ylim([-.2 .4]*4)

%legend({'PRT', 'Placebo', 'Usual Care'}, 'FontSize', 20)
legend off

% add error bars
for i=1:3
    vals = tabl.X1(tabl.group==i);
    h=errorbar(mean(vals), .8 + i/4, ste(vals)*1.96, 'horizontal', 'LineWidth', 6, 'CapSize',15,'Color', colors(i,:))
end

% print pdf
set(gcf, 'Position', [160        393        1068         162])
g = gcf;
g.PaperSize = [20 5]; % keep PaperPositionMode at manual. Might need to adjust PaperPositionSize though, or, this variable

%saveas(gcf, fullfile(basedir, 'figures', 'document_scores_1st_dim.pdf'))

%% ranksum test, effect size PRT vs PLA
clc
[p,~, stats] = ranksum(tabl.X1(tabl.group==1), tabl.X1(tabl.group==2))

g = mes(tabl.X1(tabl.group==1), tabl.X1(tabl.group==2), 'hedgesg');
g.hedgesg, g.hedgesgCi

%% ranksum test, effect size PRT vs UC
clc
[p,h, stats] = ranksum(tabl.X1(tabl.group==1), tabl.X1(tabl.group==3))

g = mes(tabl.X1(tabl.group==1), tabl.X1(tabl.group==3), 'hedgesg');
g.hedgesg, g.hedgesgCi

%% does semantic location relate to pain?
tabl.delta_pain_z = zscore(tabl.pain_avg_t2_arm_1 - tabl.pain_avg_baseline);
tabl.X1_z = zscore(tabl.X1);
fitglm(tabl, 'delta_pain_z ~ group*X1_z + age + gender', 'CategoricalVars','group')