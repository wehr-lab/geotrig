function [cricket_present]=get_cricket_present_frames(metadata, localframe, num_geoframes, filename) 

% usage: [cricket_present]=get_cricket_present_frames(metadata, localframe, num_geoframes, filename) 
%
% determine cricket_present frames, i.e. 1 in between cricket drop and
% capture frame, and 0 everywhere else

tic
fprintf('\ngetting cricket-present frames... ')
cricket_present=zeros(num_geoframes, 1);
for i=1:height(metadata)
    cricketdrop = metadata{i, 'cricketdrop'};
    captureframe = metadata{i, 'captureframe'};
    fname=metadata{i, 'filename'};
    trialframes=find(contains(filename, fname));
    cricket_present(trialframes) = (localframe(trialframes) > cricketdrop) & (localframe(trialframes) < captureframe)  ;
    %frames that are after cricket drop and before captureframe and are on this trial
   
end
fprintf(' done (%.0f sec)', toc)
