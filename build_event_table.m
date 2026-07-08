%% =========================================================================
%  BUILD_EVENT_TABLE
%  Build event_tbl from separate event-frame lists (one vector per event
%  type), looking up TrialID/AnimalID for each event via trial_id and
%  animal_id (from derive_trial_and_animal_id.m).
%
%  Assumes: contact_gain_event_frames, contact_loss_event_frames,
%  cricket_jump_event_frames, failed_approach_event_frames,
%  intercept_event_frames, rangemin_event_frames are all frame-index
%  vectors indexing into the SAME frame timeline as trial_id/animal_id.
% =========================================================================

eventFrameLists = struct( ...
    'contact_gain',    contact_gain_event_frames(:), ...
    'contact_loss',    contact_loss_event_frames(:), ...
    'cricket_jump',    cricket_jump_event_frames(:), ...
    'failed_approach', failed_approach_event_frames(:), ...
    'intercept',       intercept_event_frames(:), ...
    'range_min',       rangemin_event_frames(:));

eventTypes = fieldnames(eventFrameLists);
nTotal = sum(structfun(@numel, eventFrameLists));

EventFrame = nan(nTotal, 1);
EventType  = strings(nTotal, 1);

row = 0;
for k = 1:numel(eventTypes)
    frames = eventFrameLists.(eventTypes{k});
    n = numel(frames);
    EventFrame(row+1 : row+n) = frames;
    EventType(row+1 : row+n)  = eventTypes{k};
    row = row + n;
end

%% 1. Validate frame indices are within range --------------------------------
nFramesTotal = numel(trial_id);
badIdx = isnan(EventFrame) | EventFrame < 1 | EventFrame > nFramesTotal;
if any(badIdx)
    warning('%d / %d event(s) have out-of-range or NaN frame indices -- dropping.', ...
        nnz(badIdx), nTotal);
    EventFrame(badIdx) = [];
    EventType(badIdx)  = [];
end

%% 2. Look up TrialID / AnimalID at each event frame -------------------------
TrialID  = trial_id(EventFrame);
AnimalID = animal_id(EventFrame);

%% 3. Assemble table -----------------------------------------------------------
event_tbl = table(EventFrame, EventType, TrialID, AnimalID);
event_tbl = sortrows(event_tbl, 'EventFrame');

%% 4. Sanity checks --------------------------------------------------------------
fprintf('Built event_tbl with %d events across %d event types:\n', ...
    height(event_tbl), numel(eventTypes));
disp(groupsummary(event_tbl, 'EventType'));

% Duplicate (frame, type) rows usually indicate an upstream bug (same
% event logged twice, or two detectors firing on the same frame).
[~, uRows] = unique(event_tbl(:, {'EventFrame','EventType'}), 'rows');
nDup = height(event_tbl) - numel(uRows);
if nDup > 0
    warning('%d duplicate (EventFrame, EventType) row(s) detected -- inspect upstream detection.', nDup);
end

% Flag any event whose frame's AnimalID is NaN (e.g. unparsed filename)
nMissingAnimal = nnz(isnan(event_tbl.AnimalID));
if nMissingAnimal > 0
    warning('%d event(s) have missing AnimalID -- check derive_trial_and_animal_id.m output.', ...
        nMissingAnimal);
end
