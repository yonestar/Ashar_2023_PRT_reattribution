basedir = '/Users/yoni/Repositories/OLP4CBP/scripts/analyses/semantic_analyses';
figdir = fullfile(basedir, 'figures');
attr_wide = readtable(fullfile(basedir, 'manual_codings', 'IPQlastitem_final MERGED WIDE FORMAT.csv'));
attr_long = readtable(fullfile(basedir, 'manual_codings', 'IPQlastitem_final MERGED LONG FORMAT.csv'));

% code for time, in long form
a=categorical(attr_long.redcap_event_name);
attr_long.time = grp2idx(a);

ngram_len = 1; %1:2

%% get words from T1 into variables

attrT1 = attr_long(strcmp('t1_arm_1', attr_long.redcap_event_name),:);
attrT2 = attr_long(strcmp('t2_arm_1', attr_long.redcap_event_name),:);

% get only data from completers
attrT1_completers = attrT1(ismember(attrT1.id, attrT2.id),:);


%% get words from PRT group into variables -- completers only
PRT_T1_words = cat(1, attrT1_completers.attributed_cause_1, attrT1_completers.attributed_cause_2, attrT1_completers.attributed_cause_3);
PRT_T2_words = cat(1, attrT2.attributed_cause_1, attrT2.attributed_cause_2, attrT2.attributed_cause_3);

%% top words in PRT group at T1 and T2 -- set up for "biggest movers"

n = 10000; % to get them ALL
clc

% preproc
%preproced_T1 = preprocessTextData(PRT_T1_words);
%preproced_T2 = preprocessTextData(PRT_T2_words);

t1_ngrams = bagOfNgrams(preprocessTextData(PRT_T1_words), 'NgramLengths',ngram_len);
topkT1 = topkngrams(t1_ngrams, n);

t2_ngrams = bagOfNgrams(preprocessTextData(PRT_T2_words), 'NgramLengths',ngram_len);
topkT2 = topkngrams(t2_ngrams, n);

% for display
n=20;
topkngrams(t1_ngrams, n)
topkngrams(t2_ngrams, n)

%% words with largest increases in frequency, T1 to T2
clc
attr_changes = outerjoin(topkT1, topkT2, 'Keys', 'Ngram'); 
attr_changes.Count_topkT2(isnan(attr_changes.Count_topkT2)) = 0; % IF WORD ONLY EXISTS AT ONE TIMEPOINT, make it 0 at the other timepoint
attr_changes.Count_topkT1(isnan(attr_changes.Count_topkT1)) = 0; % IF WORD ONLY EXISTS AT ONE TIMEPOINT, make it 0 at the other timepoint

attr_changes.delta_count = attr_changes.Count_topkT2 - attr_changes.Count_topkT1; 

n=22; whcol = [4 2 5 7];
attr_changes = sortrows(attr_changes, 'delta_count', 'descend');
topk = head(attr_changes(:,whcol),n)

% find examples
topk.Ngram = topk.Ngram_topkT2;
%findSomeExamples(topk, cat(1, PRT_T1_words, PRT_T2_words))

%% words with largest decreases in frequency, T1 to T2
clc
attr_changes = sortrows(attr_changes, 'delta_count', 'ascend');
n=25; whcol = [1 2 5 7];
topk = head(attr_changes(:,whcol),n)

%% find examples
topk.Ngram = topk.Ngram_topkT1;
findSomeExamples(topk, cat(1, PRT_T1_words, PRT_T2_words))


% attr_changes = innerjoin(countsT1, countsT2, 'Keys', 'Word');
% attr_changes.delta_count = attr_changes.Count_countsT2 - attr_changes.Count_countsT1; 
% 
% attr_changes = sortrows(attr_changes, 'delta_count', 'descend');
% head(attr_changes(:,[1 4]),n)
% 
% attr_changes = sortrows(attr_changes, 'delta_count', 'ascend');
% head(attr_changes(:,[1 4]),n)

%% SUBFUNCTIONS

% important note -- this won't find ALL the examples, b/c of lemmatization
% process
function findSomeExamples(topk, words)

for i=1:height(topk)
    % take one or two words, as needed
    if topk.Ngram(i,2)==""
        myword = topk.Ngram(i,1);
    else
        myword = convertStringsToChars([topk.Ngram{i,1} ' ' topk.Ngram{i,2}]);
    end
    
    fprintf('WORD: %s\n', myword);
    words(contains(words, myword,'IgnoreCase',true))
end

end