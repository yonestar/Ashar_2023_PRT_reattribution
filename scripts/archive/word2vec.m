
%emb = fastTextWordEmbedding;
T1wordstable = readtable('T1wordstable.csv')
3columnwords = T1wordstable(:,3:5)
%words = reshape(3columnwords,[],1)


%% 

wordvectors = word2vec(emb,words)
XY = tsne(wordvectors, 'NumDimensions', 3)
numWords = numel(wordvectors);
colorData = rand(numWords,3);
figure
textscatter3(XY,wordvectors,'ColorData',colorData)
title("Word Embedding t-SNE Plot")



%% 



% to get the mean of the 3 words
all3words = cat(3, word2vec(emb, T1WordsTable.st_most_important), word2vec(emb, T1WordsTable.nd_most_important));

mean(all3words,3)

M = word2vec(emb,words)

%% clean stop words from tokenized document
newDocuments = removeStopWords(documents)



%% or try glove -- must download file first

filename = "glove.6B.300d";
if exist(filename + '.mat', 'file') ~= 2
    emb = readWordEmbedding(filename + '.txt');
    save(filename + '.mat', 'emb', '-v7.3');
else
    load(filename + '.mat')
end
v_king = word2vec(emb,'king')';
whos v_king


%% load data
load('T1WordsTable.mat');

% to get the mean of the 3 words
all3words = cat(3, word2vec(emb, T1WordsTable.st_most_important), word2vec(emb, T1WordsTable.nd_most_important));

mean(all3words,3)


%%
words= ["spine";"injury";"scoliosis";"Accident"; "emotion"; "brain"; "mind"; "stress"; "pregnancy"]
%for sub = 1:length(words)   
M = word2vec(emb,words)
XY = tsne(M, 'NumDimensions', 3)
numWords = numel(words);
colorData = rand(numWords,3);
figure
textscatter3(XY,words,'ColorData',colorData)
title("Word Embedding t-SNE Plot")

%{
 questions:

how to deal with multiple words (phrases) -- see their preprocessing stuff
is fastTextEmbedding a good pre-trained space?
is tsne appropriate, vs. pca
whether plot in the space in general vs. plot distance from predtermined
points of theoretical interest


%}