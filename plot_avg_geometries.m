function plot_avg_geometries(event_frames, event_type, mouse_spd, range, az, cricket_spd, close_fig, outputrootdir)
%usage: plot_avg_geometries(event_frames, mouse_spd, range, az, cricket_spd, close_fig, outputrootdir)
% plots event-triggered average geometries (az, range, speeds) based on events in event_frames
% appends to pdf file event-avg-geo.pdf in outputrootdir
% set close_fig to 1 to close figures immediately after exporting to pdf

winsize_geo=1*200;

win=(-winsize_geo):(+winsize_geo);
%initialize
all_az=nan(length(event_frames), length(win));
all_range=all_az;
all_mouse_spd=all_az;
all_cricket_spd=all_az;

nr=nanstd(range);
na=nanstd(az);
ns=nanstd(mouse_spd);
ncs=nanstd(cricket_spd);

p=0;offset=0;

for event_frame=event_frames(:)'; %at 200 fps
    p=p+1;
    win=(event_frame-winsize_geo):(event_frame+winsize_geo);

    all_az(p, :)=az(win)/na;
    all_range(p, :)=range(win)/nr;
    all_mouse_spd(p, :)=mouse_spd(win)/ns;
    all_cricket_spd(p, :)=cricket_spd(win)/ncs;
end

mean_az=nanmean(all_az);
mean_range=nanmean(all_range);
mean_mouse_spd=nanmean(all_mouse_spd);
mean_cricket_spd=nanmean(all_cricket_spd);
t= -winsize_geo:+winsize_geo;

mkdir(outputrootdir)
pdffilename=fullfile(outputrootdir, 'event-avg-geo.pdf');

figure
hold on
plot(t, mean_az,  'linew', 2)
plot(t, mean_range,  'linew', 2)
plot(t, mean_mouse_spd,  'linew', 2)
plot(t, mean_cricket_spd,  'linew', 2)
line([0 0], ylim)
legend('az', 'range', 'mouse speed', 'cricket speed', 'Location','northeastoutside')
title(sprintf('mean geometries aligned to %s', event_type), 'interpreter', 'none')
set(gcf, 'pos' , [1000         818         819         420])
xlabel(sprintf('frames, n=%d events', length(event_frames)))
ylabel('z normalized')
exportgraphics(gcf,pdffilename,'Append',true)
if close_fig
    close
end
