function [lightdark light dark] = get_lightdark(filename, localframe, num_geoframes, metadata);

% usage: [lightdark light dark] = get_lightdark(filename, localframe, num_geoframes, metadata);
% determine light or dark condition for every frame
% lightdark is 1 for light, 2 for dark, 0 for other
% light and dark are logicals (1 when true)

tic
fprintf('\ngetting light/dark frames... ')
lightdark=zeros(num_geoframes, 1);
for i=1:height(metadata)
    condition = metadata{i, 'condition'};
    fname=metadata{i, 'filename'};
    trialframes=find(contains(filename, fname));
    lightdark(trialframes) = condition ;
   
end
fprintf(' done (%.0f sec)', toc)
light=lightdark==1;
dark=lightdark==2;
