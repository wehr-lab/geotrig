function capture_frames=detect_capture(metadata, filename);

% usage: 
%returns capture frame (the one at the end of each trial) in absolute frames
capture_frames=[];
fprintf('detecting capture frames... ')

nbytes =  fprintf('%.1f%%', 0);
for i=1:height(metadata)
   fprintf(repmat('\b',1,nbytes)); nbytes =  fprintf('%.1f%%', 100*i/height(metadata)); 
    fn=metadata{i, 'filename'};
    capframe=metadata{i, 'captureframe'};
    startframe=(find(contains(filename, fn), 1));
    if ~isempty(startframe)
        capture_frames=[capture_frames capframe+startframe];
    end
end
capture_frames=capture_frames(:);            
fprintf('done. found %d capture frames\n', length(capture_frames))

