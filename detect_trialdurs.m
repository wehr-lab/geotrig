function trial_durations=detect_trialdurs(metadata, filename);

% usage:  trial_durations=detect_trialdurs(metadata, filename);
%
%returns trial_durations (num frames between cricket drop and captureframe)

capture_frames=[];
fprintf('detecting trial_durations... ')

nbytes =  fprintf('%.1f%%', 0);
trial_durations=[];
for i=1:height(metadata)
    fprintf(repmat('\b',1,nbytes)); nbytes =  fprintf('%.1f%%', 100*i/height(metadata));
    fn=metadata{i, 'filename'};
    dropframe=metadata{i, 'cricketdrop'};
    capframe=metadata{i, 'captureframe'};
    if capframe~=0 & dropframe~=0
        trial_durations=[trial_durations capframe-dropframe];
    end
end
trial_durations=trial_durations(:);
fprintf('done. found %d trial_durations (min %d, max %d, mean %.1f)\n', length(trial_durations), min(trial_durations), max(trial_durations), nanmean(trial_durations))

