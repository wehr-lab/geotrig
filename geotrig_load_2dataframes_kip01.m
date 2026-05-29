function [az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata] = geotrig_load_dataframe_kip01(datapath, dataframeFN, metadataFN);

% usage: [az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata] = geotrig_load_dataframe(datapath, dataframeFN);
% dataframeFN: name of monolithic file with all mice, all trials, all frames. For example 'geos_allmice_alltrials_wfames.csv'; 
% loads the entire geometries dataframe
% it's slow, could take up to a few minutes
%
% example call: 
% datapath= '/Volumes/Projects/PreyCapture/ZI_A1Activation';
% dataframeFN = 'geos_allmice_alltrials_wfames.csv'; 
% [az, range, localframe, cricket_spd, mouse_spd, localframe, filename, num_geoframes, metadata] = geotrig_load_dataframe(datapath, dataframeFN);

cd(datapath)

tic
if 0
    fprintf('detecting import opts ...\n')
    opts = detectImportOptions(dataframeFN, 'ExpectedNumVariables', 6, 'ReadVariableNames',true, 'ReadRowNames',false, 'NumHeaderLines',1, 'Delimiter',',', ...
        'EmptyColumnRule', 'read', 'ExtraColumnsRule', 'addvars',  'VariableNamesLine', 1);
    opts.DataLines=[2 inf];
    %    opts.VariableTypes={'char','double','double','double','double','double'};
    %this seems to break when switching to different dataframe csvs, it seems
    %to be fooled by nan columns, not sure why overriding opts doesn't always
    %work (because it only looke at the first 200 lines!)
    % columns: session#, filename, localframe, dist, cricket_spd,	mouse_spd, az
else
    fprintf('delimiting import opts ...\n')
    opts =  delimitedTextImportOptions('VariableNames', {'filename', 'localframe', 'dist', 'cricket_spd', 'mouse_spd', 'az'} ,...
        'VariableTypes',{'char','double','double','double','double','double'}, 'DataLines', [2] ,'VariableNamesLine', [1],...
        'Delimiter',',');
end
elapsed_time=toc;


fprintf('done %.0fs\n', elapsed_time)

fprintf('reading entire dataframe...\n')
tic
dataframe = readtable(dataframeFN,opts, ReadRowNames=false, ReadVariableNames=true);

metadata=readtable(metadataFN,'Delimiter',',','ReadVariableNames',true, 'ReadRowNames',false);

num_geoframes=height(dataframe);
az=dataframe.az;
range=dataframe.dist;
%localframe=dataframe.localframe;
cricket_spd=dataframe.cricket_spd;
mouse_spd=dataframe.mouse_spd;
localframe=dataframe.localframe;
indStarts = find(localframe==0);    % these are the row indices into the dataframe for the start of each test
filename=dataframe.filename;
elapsed_time=toc;
fprintf('done. read %d frames in %.0fs\n', num_geoframes, elapsed_time)

