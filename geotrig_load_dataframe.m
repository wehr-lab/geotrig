function [az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata] = geotrig_load_dataframe(datapath, dataframefilename);

% usage: [az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata] = geotrig_load_dataframe(datapath, dataframefilename);
% dataframefilename: name of monolithic file with all mice, all trials, all frames. For example 'geos_allmice_alltrials_wfames.csv'; 
% loads the entire geometries dataframe
% it's slow, could take up to a few minutes
%
% example call: 
% datapath= '/Volumes/Projects/PreyCapture/ZI_A1Activation';
% dataframefilename = 'geos_allmice_alltrials_wfames.csv'; 
% [az, range, localframe, cricket_spd, mouse_spd, localframe, filename, num_geoframes, metadata] = geotrig_load_dataframe(datapath, dataframefilename);

fprintf('\nstart time %s \n detecting import opts ...', datestr(now))
cd(datapath)

% % it requires a bunch of opts specifications to read correctly
% opts = detectImportOptions(dataframefilename, 'ExpectedNumVariables', 6, 'ReadVariableNames',true, 'ReadRowNames',false, 'NumHeaderLines',1, 'Delimiter',',', ...
%     'ExtraColumnsRule', 'addvars',  'VariableNamesLine', 1);
% opts.DataLines=[2 inf];
% %opts.VariableNames{7}='az';
% opts.VariableTypes={'double','char','double','double','double','double'};       %,'double'};
% 
% %this seems to break when switching to different dataframe csvs, it seems
% %to be fooled by nan columns, not sure why overriding opts doesn't always work 
% % columns: session#, filename, localframe, dist, cricket_spd,	mouse_spd, az
%%% the problem here is that MATLAB only looks in the first 200 lines of
% the file to make its determination of #columns...

% This works if you know what the file structure is.
fprintf('delimiting import opts ...\n')
opts =  delimitedTextImportOptions('VariableNames', ...
    {'sessionID', 'filename', 'localframe', 'dist', 'cricket_spd', 'mouse_spd', 'az'} ,...
    'VariableTypes',{'char','char','double','double','double','double','double'}, ...
    'DataLines', [2] ,'VariableNamesLine', [1],'Delimiter',',');
fprintf('done')

fprintf('\nreading entire dataframe...')
tic
dataframe = readtable(dataframefilename,opts, ReadRowNames=true, ReadVariableNames=true);

metadata=readtable('metadata_alltrials.csv','Delimiter',',');

num_geoframes=height(dataframe);
az=dataframe.az;
range=dataframe.dist;
cricket_spd=dataframe.cricket_spd;
mouse_spd=dataframe.mouse_spd;
localframe=dataframe.localframe;
filename=dataframe.filename;
elapsed_time=toc;
fprintf('done. read %d frames in %.1fmin', num_geoframes, elapsed_time/60)
%keyboard

%%%%%%%%%%%%%%%

% 6.1.26
% opts used to be like this, not sure why it stopped working
% 
% opts =  delimitedTextImportOptions('VariableNames', ...
%     {'filename', 'localframe', 'dist', 'cricket_spd', 'mouse_spd', 'az'} ,...
%     'VariableTypes',{'char','double','double','double','double','double'}, ...
%     'DataLines', [2] ,'VariableNamesLine', [1],'Delimiter',',');
% fprintf('done')

% % for troubleshooting, the dataframe is so huge that it's convenient to use a small temp version
% temp_dataframefilename='temp_geos_allmice.csv'
% fid=fopen(dataframefilename, 'r')
% fid2=fopen('temp_geos_allmice.csv', 'w')
% for i=1:1000
%     l=fgetl(fid)
%     fprintf(fid2, '%s\n', l)
%     l=fgetl(fid)
%     fprintf(fid2, '%s', l)
% end
% fclose(fid)
% fclose(fid2)




