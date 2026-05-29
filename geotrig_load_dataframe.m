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

fprintf('\ndetecting import opts ...')
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
opts =  delimitedTextImportOptions('VariableNames', {'filename', 'localframe', 'dist', 'cricket_spd', 'mouse_spd', 'az'} ,...
    'VariableTypes',{'char','double','double','double','double','double'}, 'DataLines', [2] ,'VariableNamesLine', [1],...
    'Delimiter',',');
fprintf('done')

fprintf('\nreading entire dataframe...')
tic
dataframe = readtable(dataframefilename,opts, ReadRowNames=true, ReadVariableNames=true);

metadata=readtable('metadata_alltrials.csv','Delimiter',',');

num_geoframes=height(dataframe);
az=dataframe.az;
range=dataframe.dist;
%localframe=dataframe.localframe;
cricket_spd=dataframe.cricket_spd;
mouse_spd=dataframe.mouse_spd;
localframe=dataframe.localframe;
filename=dataframe.filename;
elapsed_time=toc;
fprintf('done. read %d frames in %.0fs', num_geoframes, elapsed_time)

