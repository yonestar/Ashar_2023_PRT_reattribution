basedir = '/Users/yoni/Repositories/OLP4CBP/scripts/analyses/semantic_analyses';

f1 = readtable(fullfile(basedir, 'IPQlastitem_final_deid DK.csv'));
f2 = readtable(fullfile(basedir, 'IPQlastitem_final_deid JK.csv'));

clc

% drop the stars before reconciling
f1 = dropstars(f1);
f2 = dropstars(f2);

%% how many total diffs?
diffs = setdiff(f1,f2);

%% reconcile
numdiffs = 0;
for i=1:height(f1)
    
    if ~isequal(f1{i, end-2:end}, f2{i, end-2:end})
        [f1(i,:); f2(i,:)]
        numdiffs = numdiffs+1;
    end
end

%% subfunctions
function tabl = dropstars(tabl)

    for i=5:7
        wh = strcmp(tabl{:,i}, 'MB*'); 
        tabl(wh, i) = {'MB'};
        
        wh = strcmp(tabl{:,i}, 'STR*'); 
        tabl(wh, i) = {'STR'};
    end
end