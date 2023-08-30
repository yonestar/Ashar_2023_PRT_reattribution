basedir = '/Users/yoni/Repositories/Ashar_2023_PRT_reattribution/data/manual_codings';

f1 = readtable(fullfile(basedir, 'IPQlastitem_final_deid categories YA.csv'));
f2 = readtable(fullfile(basedir, 'IPQlastitem_final_deid categories EL.csv'));

clc

% convert all to lower
f1.cat1 = lower(f1.cat1);
f1.cat2 = lower(f1.cat2);
f1.cat3 = lower(f1.cat3);

f2.cat1 = lower(f2.cat1);
f2.cat2 = lower(f2.cat2);
f2.cat3 = lower(f2.cat3);


%% how many total diffs?
clc
numdiffs = 0;
for i=1:height(f1)%  -- I reconciled up through 150. Notes to self: MB vs STR ratings are highly reliable between coders, but categorizations are quite variable and/or require much more training and clear boundaries and procedures when attributions are multi-modal/complex
    
    if ~isequal(f1{i, end-2:end}, f2{i, end-2:end})
        disp([f1(i,:); f2(i,:)])
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