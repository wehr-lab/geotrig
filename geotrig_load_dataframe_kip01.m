function [az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata, dataframe] = geotrig_load_dataframe_kip01(datapath, dataframeFN1, metadataFN1, dataframeFN2, metadataFN2);

% usage: [az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata, dataframe] = geotrig_load_dataframe(datapath, dataframeFN);
% dataframeFN: name of monolithic file with all mice, all trials, all frames. For example 'geos_allmice_alltrials_wfames.csv'; 
% loads the entire geometries dataframe
% it's slow, could take up to a few minutes
%
% example call: 
% datapath= '/Volumes/Projects/PreyCapture/ZI_A1Activation';
% dataframeFN = 'geos_allmice_alltrials_wfames.csv'; 
% [az, range, localframe, cricket_spd, mouse_spd, localframe, filename, num_geoframes, metadata] = geotrig_load_dataframe(datapath, dataframeFN);

cd(datapath)

fprintf('delimiting import opts ...\n')
opts =  delimitedTextImportOptions('VariableNames', {'filename', 'localframe', 'dist', 'cricket_spd', 'mouse_spd', 'az'} ,...
    'VariableTypes',{'char','double','double','double','double','double'}, 'DataLines', [2] ,'VariableNamesLine', [1],...
    'Delimiter',',');

fprintf('reading entire dataframe1...\n')
tic
dataframe = readtable(dataframeFN1,opts, ReadRowNames=false, ReadVariableNames=true);

metadata = readtable(metadataFN1,'Delimiter',',','ReadVariableNames',true, 'ReadRowNames',false);
num_geoframes=height(dataframe);
elapsed_time=toc;

fprintf('done. read %d frames in %.0fs\n', num_geoframes, elapsed_time)

if nargin ==5
    fprintf('reading entire dataframe2...\n')
    tic
    dataframe2 = readtable(dataframeFN2,opts, ReadRowNames=false, ReadVariableNames=true);

    metadata2 = readtable(metadataFN2,'Delimiter',',','ReadVariableNames',true, 'ReadRowNames',false);
    num_geoframes=height(dataframe2);
    elapsed_time=toc;

    fprintf('done. read %d frames in %.0fs\n', num_geoframes, elapsed_time)

    % vertically concatenate dataframe and metadata
    dataframe = [dataframe; dataframe2];
    metadata = [metadata; metadata2];
end

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


