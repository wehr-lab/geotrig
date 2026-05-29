function [laseron] = get_laseron(filename, localframe, num_geoframes, metadata);

% usage: [ laseron] = get_laser_on(filename, localframe, num_geoframes, metadata);
% determine laser on or off condition for every frame
% laseron is 1 for on, 0 for off, it's a logical 

tic
fprintf('\ngetting laseron frames... ')
laseron=zeros(num_geoframes, 1);
for i=1:height(metadata)
    laser_value = metadata{i, 'laser_value'};
    fname=metadata{i, 'filename'};
    trialframes=find(contains(filename, fname));
    laseron(trialframes) = laser_value ;
   
end
fprintf(' done (%.0f sec)', toc)
laseron(isnan(laseron))=0;
laseron=logical(laseron);