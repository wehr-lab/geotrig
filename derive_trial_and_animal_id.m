%% =========================================================================
%  DERIVE_TRIAL_AND_ANIMAL_ID
%  Extract trial_id and animal_id from `filename` (numframes x1 cell array
%  of strings), where each unique filename = one trial, and each filename
%  contains a substring like "mouse-0897" identifying the animal.
% =========================================================================

filename = string(filename(:));   % ensure Nx1 string array (accepts cellstr too)

%% 1. trial_id: integer per unique filename (each filename = one trial) ----
[uniqueFiles, ~, trial_id] = unique(filename, 'stable');
trial_id = trial_id(:);

%% 2. animal_id: parse "mouse-####" substring, use the NUMBER as the ID ----
% Using the actual mouse number (rather than an order-of-appearance index)
% keeps animal_id stable and meaningful across separate scripts/sessions.
tok = regexpi(filename, '(?<=mouse-)\d+', 'match', 'once');
missing = (tok == "");
if any(missing)
    warning('%d / %d frames had no "mouse-####" match in filename.', ...
        nnz(missing), numel(filename));
end
animal_id = nan(numel(filename), 1);
animal_id(~missing) = str2double(tok(~missing));

%% 3. Sanity checks ----------------------------------------------------------
% Each trial (= each unique filename) should map to exactly one animal.
% If not, the filename parsing or trial/animal assumption needs a look.
g = findgroups(trial_id);
nAnimalsPerTrial = splitapply(@(x) numel(unique(x(~isnan(x)))), animal_id, g);
if any(nAnimalsPerTrial > 1)
    badTrials = uniqueFiles(nAnimalsPerTrial > 1);
    warning('%d trial(s) contain more than one animal_id. Offending filenames:', numel(badTrials));
    disp(badTrials);
end

fprintf('Parsed %d frames -> %d unique trials, %d unique animals.\n', ...
    numel(filename), numel(uniqueFiles), numel(unique(animal_id(~isnan(animal_id)))));
