function GenerateEventVideoClips(event_frames, event_type, filename, localframe, allnumframes, max_num_clips, outputrootdir)
% generates a short video clip for each event in event_frames, saves to an 'event' folder inside of outputrootdir
% usage: GenerateEventVideoClips(event_frames, event_type, filename, localframe, allnumframes, max_num_clips, outputrootdir)
% event_frames: align clips to these frames
% event_type: name for folder
% filename, localframe: vars from dataframe 
% allnumframes: use a single value (like 400) to use the same numframes for all events (fixed duration)
% if each event has a different duration, allnumframes should be the durations for all events in event_frames
% max_num_clips: only generate this many clips. Should be an integer multiple of 16 for tiling into a 4x4 mosaic

outputdir= fullfile(outputrootdir,event_type);
outputdir= erase(outputdir,whitespacePattern)   ;

mkdir(outputdir)
cd(outputdir)
delete('*.mp4')
if length(event_frames)>max_num_clips
    event_frames=event_frames(1:max_num_clips);
end

tic
for j=1:length(event_frames)
    event_frame=event_frames(j); %at 200 fps

    fprintf('\nevent_type %s, event_frame %d, localframe %d,', event_type, event_frame, localframe(event_frame))
    % fprintf(' localframe %d, filename %s', localframe(event_frame), filename{event_frame})
    movdir=dir(filename{event_frame}).folder;
    cd(movdir)
    movfiledir=dir('Sky*labeled.mp4');
    movfilename=fullfile(movfiledir.folder, movfiledir.name); %abs path
    outputfilename=sprintf('clip_%d', j);
    try
        if length(allnumframes)>1 %each event has a different duration
            numframes=allnumframes(j);
            keyframe=-1;
            startframe=localframe(event_frame);
        else %use the same numframes for all events (fixed duration)
            numframes=allnumframes;
            startframe=localframe(event_frame)-round(numframes/2);
            if startframe<1 startframe=1;end
            keyframe=round(numframes/2);
        end
        ExtractVideoClip_ffmpeg(movfilename, startframe, numframes, outputdir, outputfilename, keyframe)
        fprintf('\n%d/%d wrote %s (%d frames, %.0f sec)', j, length(event_frames), outputfilename,numframes, toc)
    catch
        fprintf('\n error from ExtractVideoClip %d %s %s', j,  outputfilename, lasterr)
    end
end
fprintf('\n\nwrote %d clips in %.0f sec', length(event_frames), toc)

cd(outputdir)
TileVidClips_ffmpeg


% %% TileVidClips for events
% for e=1:length(event_type)
%     fprintf('\n%d/%d %s', e, length(event_type), event_type{e});
%     %j=0;
%     outputdir= sprintf('/Volumes/Projects/PreyCapture/ZI_A1Activation/geo-trig-analysis-output/%s',event_type{e});
%     outputdir= erase(outputdir,whitespacePattern)   ;
%     %outputdir=[outputdir, '1']
%     %mkdir(outputdir)
%     cd(outputdir)
%     TileVidClips_ffmpeg
% end