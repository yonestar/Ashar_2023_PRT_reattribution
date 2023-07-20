% 1st varargin = keep empty documents?

function documents = preprocessTextData(textData, varargin)

if nargin > 1
    keepEmpty = true;
else
    keepEmpty = false; 
end

% erratic behavior when text is mix of upper and lower case. not sure why.
% "Neural Pathways" and "neural pathways" yields different results.
% Capitalized words may be processed differenltly?
textData = lower(textData);

% Tokenize the text.
documents = tokenizedDocument(textData);

% spell check
documents = correctSpelling(documents);

% Lemmatize the words. To improve lemmatization, first use 
% addPartOfSpeechDetails.
documents = addPartOfSpeechDetails(documents);
documents = normalizeWords(documents,'Style','lemma');

% Erase punctuation.
documents = erasePunctuation(documents);

% Remove a list of stop words.
documents = removeStopWords(documents);

% Remove words with 2 or fewer characters, and words with 15 or more
% characters.
documents = removeShortWords(documents,2);
documents = removeLongWords(documents,20);

if ~keepEmpty, documents = removeEmptyDocuments(documents); end

end