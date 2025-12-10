function [target_loss_event_frames]=detect_target_loss(range, az, contact, metadata, localframe, num_geoframes, filename) 

% usage: target_loss_event_frames=detect_target_loss(range, az, contact, metadata, localframe, num_geoframes, filename); 
% detects target_loss events, defined by loss of contact followed by a
% window during which az and range exceed thresholds


% params:
target_loss_range_thresh=20; %range definition of target_loss, cm
target_loss_az_thresh=45; %az definition of target_loss, deg
winsize=2*200; %in frames (seconds*200fps)

tic
fprintf('\n')
nbytes =  fprintf('%.1f%%', 0);
target_loss_event_frames=[];
for f=1:num_geoframes-winsize-1
    if ~mod(f, 100000), fprintf(repmat('\b',1,nbytes)); nbytes =  fprintf('%.1f%%', 100*f/num_geoframes); end

    if contact(f)==1 & contact(f+1)==0 %filter by loss of contact
        cricketdrop = metadata{contains(metadata.filename, filename{f}), 'cricketdrop'};
        captureframe = metadata{contains(metadata.filename, filename{f}), 'captureframe'};
        if localframe(f)>cricketdrop & localframe(f)<captureframe         %exclude frames before cricket drop or after capture
            if median(range(f+1:f+1+winsize))>=target_loss_range_thresh ...
                    & nanmedian(az(f+1:f+1+winsize))>=target_loss_az_thresh
                target_loss_event_frames=[target_loss_event_frames f];
            end
        end
    end
end
fprintf(repmat('\b',1,nbytes));
fprintf('%.1f%%', 100*f/num_geoframes);
fprintf(' done (%.0f sec)', toc)


