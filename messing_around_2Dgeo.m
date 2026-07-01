%messing around with visualizing 2D & 3D geometry distributions and search/pursuit

%% load geos from dataframe
if 0
    dataroot= '/Volumes/Projects/PreyCapture/ZIActivation';
    dataframefilename = 'geos_allmice_alltrials_wfames.csv';
    outputrootdir= '/Volumes/Projects/PreyCapture/ZIActivation/geo-trig-analysis-output';
    [az, range, localframe, cricket_spd, mouse_spd, filename, num_geoframes, metadata] = geotrig_load_dataframe(dataroot, dataframefilename);
    cd '/Users/wehr/Documents/Projects/Molly Shallow'

    cricket_present=get_cricket_present_frames(metadata, localframe, num_geoframes, filename);
    [lightdark light dark] = get_lightdark(filename, localframe, num_geoframes, metadata);
    [laseron] = get_laseron(filename, localframe, num_geoframes, metadata);
end


%% detect events
cricket_jump_event_frames=detect_cricketjump(cricket_spd, metadata, localframe, filename);
rangemin_event_frames=detect_rangemin(range, metadata, localframe, filename);
[contact, contact_gain_event_frames, contact_loss_event_frames]=detect_contact(range, cricket_present);
[target_loss_event_frames]=detect_target_loss(range, az, contact, metadata, localframe, num_geoframes, filename);
[chase, chase_start_frames, chase_end_frames, chase_durs]=detect_chase(cricket_present, mouse_spd, range, az);
[pause, pause_start_frames, pause_end_frames, pause_durs]=detect_pause(cricket_present, mouse_spd);
[wander, wander_start_frames, wander_end_frames, wander_durs]=detect_wander(cricket_present, mouse_spd, range, az);
[stalk, stalk_start_frames, stalk_end_frames, stalk_durs]=detect_stalk(cricket_present, mouse_spd, cricket_spd, range, az);
[approach, approach_start_frames, approach_end_frames, approach_durs, first_approach_frames]=detect_approach(cricket_present, mouse_spd, az);



%thigmo is not in the dataframe in this dataroot

Azimuthedges=linspace(0, 180, 50);
% Rangeedges=linspace(0, 25, 50);
Rangeedges=linspace(0, 60, 50);
Speededges=linspace(0, 80, 50);

tiledlayout(1,3, "TileSpacing","compact")
nexttile
histmat=histcounts2(az, range, Azimuthedges, Rangeedges);
imagesc(Azimuthedges,Rangeedges,sqrt(histmat'));
shading interp
set(gca, 'ydir', 'norm')
xlabel('Azimuth, degrees')
ylabel('Range, cm')
title('all frames')
set(gcf, "Position", [440 880 1300 420 ])

nexttile
histmat=histcounts2(mouse_spd, range, Speededges, Rangeedges);
imagesc(Rangeedges,Speededges,sqrt(histmat));
shading interp
set(gca, 'ydir', 'norm')
xlabel('Range, cm')
ylabel('Speed, cm/s')
title('all frames')

nexttile
% histmat=histcounts2(mouse_spd, az, Speededges, Azimuthedges);
histmat=histcounts2(mouse_spd, az, Speededges, Azimuthedges);
imagesc(Azimuthedges,Speededges,sqrt(histmat));
shading interp
set(gca, 'ydir', 'norm')
xlabel('Azimuth, degrees')
ylabel('Speed, cm/s')
title('all frames')
shg

%% plot light vs dark 1-D
histnorm='probability';
azhistdark=histcounts(az(dark), Azimuthedges, 'norm', histnorm);
rangehistdark=histcounts(range(dark), Rangeedges, 'norm', histnorm);
spdhistdark=histcounts(mouse_spd(dark), Speededges, 'norm', histnorm);
azhistlight=histcounts(az(light), Azimuthedges, 'norm', histnorm);
rangehistlight=histcounts(range(light), Rangeedges, 'norm', histnorm);
spdhistlight=histcounts(mouse_spd(light), Speededges, 'norm', histnorm);
figure
tiledlayout(2,3, "TileSpacing","compact")
nexttile
plot(histx(Rangeedges), rangehistdark, histx(Rangeedges), rangehistlight)
ylabel('probability')
xlabel('Range')
box off
nexttile
plot(histx(Speededges), spdhistdark, histx(Speededges), spdhistlight)
xlabel('Mouse Speed')
box off
nexttile
plot(histx(Azimuthedges), azhistdark, histx(Azimuthedges), azhistlight)
xlabel('Azimuth')
legend('dark', 'light')
box off
xlim([0 180])


% figure
% tiledlayout(1,3, "TileSpacing","compact")
nexttile
x = histx(Rangeedges);
y =  rangehistdark- rangehistlight;
% 2. Split the data
y_pos = y; y_pos(y_pos < 0) = 0; % Keep only positive values
y_neg = y; y_neg(y_neg > 0) = 0; % Keep only negative values
hold on;
area(x, y_pos, 'FaceColor', 'r', 'EdgeColor', 'none', 'DisplayName', 'Positive');
area(x, y_neg, 'FaceColor', 'b', 'EdgeColor', 'none', 'DisplayName', 'Negative');
yline(0, 'k'); % Adds a reference line at zero
hold off;
xlabel('Range')
ylabel('dark-light')
box off

nexttile
x = histx(Speededges);
y =  spdhistdark- spdhistlight;
% 2. Split the data
y_pos = y; y_pos(y_pos < 0) = 0; % Keep only positive values
y_neg = y; y_neg(y_neg > 0) = 0; % Keep only negative values
hold on;
area(x, y_pos, 'FaceColor', 'r', 'EdgeColor', 'none', 'DisplayName', 'Positive');
area(x, y_neg, 'FaceColor', 'b', 'EdgeColor', 'none', 'DisplayName', 'Negative');
yline(0, 'k'); % Adds a reference line at zero
xlabel('Mouse Speed')
box off

nexttile
x = histx(Azimuthedges);
y =  azhistdark- azhistlight;
% 2. Split the data
y_pos = y; y_pos(y_pos < 0) = 0; % Keep only positive values
y_neg = y; y_neg(y_neg > 0) = 0; % Keep only negative values
hold on;
area(x, y_pos, 'FaceColor', 'r', 'EdgeColor', 'none', 'DisplayName', 'Positive');
area(x, y_neg, 'FaceColor', 'b', 'EdgeColor', 'none', 'DisplayName', 'Negative');
yline(0, 'k'); % Adds a reference line at zero
hold off;
xlabel('Azimuth')
box off
xlim([0 180])


%% test for independence:

histmatAR=histcounts2(az, range, Azimuthedges, Rangeedges);
histmatSR=histcounts2(mouse_spd, range, Speededges, Rangeedges);
histmatSA=histcounts2(mouse_spd, az, Speededges, Azimuthedges);
azhist=histcounts(az, Azimuthedges);
rangehist=histcounts(range, Rangeedges);
spdhist=histcounts(mouse_spd, Speededges);

% Compute Spearman's correlation coefficient and p-value (on underlying
% variables, not hists)
[rhoAR, p_valAR] = corr(az(:), range(:), 'Type', 'Spearman', 'rows', 'pairwise')
fprintf('\nAzimuth-Range: Spearman Rho: %.4f, P-value: %.4e (e%.2f)\n', rhoAR, p_valAR, log10(p_valAR));
[rhoSA, p_valSA] = corr(mouse_spd(:), az(:), 'Type', 'Spearman', 'rows', 'pairwise')
fprintf('\nSpeed-Azimuth: Spearman Rho: %.4f, P-value: %.4e (e%.2f)\n', rhoSA, p_valSA, log10(p_valSA));
[rhoSR, p_valSR] = corr(mouse_spd(:), range(:), 'Type', 'Spearman', 'rows', 'pairwise')
fprintf('\nSpeed-Range: Spearman Rho: %.4f, P-value: %.4e (e%.2f)\n', rhoSR, p_valSR, log10(p_valSR));


figure
tiledlayout(3,4, "TileSpacing","compact")

nexttile
imagesc(Rangeedges, Azimuthedges,(histmatAR));
colorbar;
title('hist');
xlabel('Range');
ylabel('Azimuth');
set(gca, 'ydir', 'norm')

% 1. Calculate total samples to normalize the counts
total_samplesAR = sum(histmatAR(:));

% 2. Calculate the expected counts (product of marginals / total)
expectedAR = (azhist(:) * rangehist(:)') / total_samplesAR;

% 3. Calculate the residuals
residualsAR = histmatAR - expectedAR;

% 4. Visualize the result
nexttile
imagesc(Rangeedges, Azimuthedges, residualsAR);
colorbar;
%title('Residual Plot (Observed - Expected)');
title(sprintf('Residual Plot (Observed - Expected)\nAzimuth-Range: Spearman Rho: %.4f, P-value: %.4e\n', rhoAR, p_valAR))
xlabel('Range');
ylabel('Azimuth');
set(gca, 'ydir', 'norm')

nexttile
imagesc(Rangeedges, Azimuthedges, expectedAR);
colorbar;
title('Expected');
xlabel('Range');
ylabel('Azimuth');
set(gca, 'ydir', 'norm')

nexttile
dependencyFactorAR = histmatAR ./ ((azhist(:) * rangehist(:)') / sum(histmatAR(:)));
imagesc(Rangeedges, Azimuthedges, dependencyFactorAR);
colorbar;
clim([-1 3]); % Adjust limits to focus on the area around 1 (independence)
set(gca, 'ydir', 'norm')
title('Dependency Factor');
xlabel('Range');
ylabel('Azimuth');

%%%%%
% Speededges, Rangeedges

nexttile
imagesc(Rangeedges, Speededges,(histmatSR));
colorbar;
title('hist');
xlabel('Range');
ylabel('Speed');
set(gca, 'ydir', 'norm')

% 1. Calculate total samples to normalize the counts
total_samplesSR = sum(histmatSR(:));

% 2. Calculate the expected counts (product of marginals / total)
expectedSR = (azhist(:) * rangehist(:)') / total_samplesSR;

% 3. Calculate the residuals
residualsSR = histmatSR - expectedSR;

% 4. Visualize the result
nexttile
imagesc(Rangeedges, Speededges, residualsSR);
colorbar;
% title('Residual Plot (Observed - Expected)');
title(sprintf('Residual Plot (Observed - Expected)\nRange-Speed: Spearman Rho: %.4f, P-value: %.4e\n', rhoSR, p_valSR))
xlabel('Range');
ylabel('Speed');
set(gca, 'ydir', 'norm')

nexttile
imagesc(Rangeedges, Speededges, expectedSR);
colorbar;
title('Expected');
xlabel('Range');
ylabel('Speed');
set(gca, 'ydir', 'norm')

nexttile
dependencyFactorSR = histmatSR ./ ((azhist(:) * rangehist(:)') / sum(histmatSR(:)));
imagesc(Rangeedges, Speededges, dependencyFactorSR);
colorbar;
clim([-1 3]); % Adjust limits to focus on the area around 1 (independence)
set(gca, 'ydir', 'norm')
title('Dependency Factor');
xlabel('Range');
ylabel('Speed');

%%% Azimuthedges, Speededges

nexttile
imagesc(Azimuthedges, Speededges,(histmatSA));
colorbar;
title('hist');
xlabel('Azimuth');
ylabel('Speed');

set(gca, 'ydir', 'norm')

% 1. Calculate total samples to normalize the counts
total_samplesSA = sum(histmatSA(:));

% 2. Calculate the expected counts (product of marginals / total)
expectedSA = (azhist(:) * rangehist(:)') / total_samplesSA;

% 3. Calculate the residuals
residualsSA = histmatSA - expectedSA;

% 4. Visualize the result
nexttile
imagesc(Azimuthedges, Speededges, residualsSA);
colorbar;
%title('Residual Plot (Observed - Expected)');
title(sprintf('Residual Plot (Observed - Expected)\nSpeed-Azimuth: Spearman Rho: %.4f, P-value: %.4e\n', rhoSA, p_valSA))
xlabel('Azimuth');
ylabel('Speed');
set(gca, 'ydir', 'norm')

nexttile
imagesc(Azimuthedges, Speededges, expectedSA);
colorbar;
title('Expected');
xlabel('Azimuth');
ylabel('Speed');
set(gca, 'ydir', 'norm')

nexttile
dependencyFactorSA = histmatSA ./ ((azhist(:) * rangehist(:)') / sum(histmatSA(:)));
imagesc(Azimuthedges, Speededges, dependencyFactorSA);
colorbar;
clim([-1 3]); % Adjust limits to focus on the area around 1 (independence)
set(gca, 'ydir', 'norm')
title('Dependency Factor');
xlabel('Azimuth');
ylabel('Speed');

set(gcf, 'pos', [720          92        1634        1240])

%%%%%%%%%%%%%%
%% get cricket moving boolean

% params:
cricket_moving_speed_thresh=1.5; %  cm/s
cricket_moving_winsize=.5*200; %in frames (seconds*200fps)

tic
fprintf('\ndetecting cricket_moving... ')

medfilt_cricket_spd = medfilt1(cricket_spd, cricket_moving_winsize, 'omitnan');
%medfilt is fast (<<1s) and by doing it here, we can use a cricket_moving-detection specific winsize
figure
[n,x]=hist(log(medfilt_cricket_spd), 100);
plot(x, n)
xlabel('log cricket speed')
line(log(cricket_moving_speed_thresh)*[1 1 ], ylim)
text(log(cricket_moving_speed_thresh), mean(ylim), sprintf('cricket-moving-speed-thresh=%.1f cm/s',cricket_moving_speed_thresh))

cricket_moving=zeros(size(cricket_spd)); %initialize to zero

% 1. Speed condition
spd_condition = medfilt_cricket_spd > cricket_moving_speed_thresh;

% 4. Combine all conditions (Boolean array)
cricket_moving = spd_condition & cricket_present;
cricket_still = ~spd_condition & cricket_present; %explicit because ~cricket_moving would include cricket not present
fprintf(' done (%.1f sec)', toc)

%%%%%%%%%%%%%%%%%%%%%%%%
%%  plot cricket moving joint probs pcolor

figure
tl=tiledlayout(4,3, "TileSpacing","compact") ;
set(gcf, "Position", [-1396.00       -958.00       1300.00       1547.00])
title(tl, 'no events or states')
i=0;

nexttile
histmat=histcounts2(az(logical(dark.*cricket_moving)), range(logical(dark.*cricket_moving)), Azimuthedges, Rangeedges);
pcolor(histx(Azimuthedges),histx(Rangeedges),sqrt(histmat'));
i=i+1;cl(i,:)=clim;
shading interp
xlabel('Azimuth, degrees')
ylabel('Range, cm')
title('cricket moving, Dark')
hold on;

nexttile
histmat=histcounts2(mouse_spd(logical(dark.*cricket_moving)), range(logical(dark.*cricket_moving)), Speededges, Rangeedges);
pcolor(histx(Rangeedges), Speededges,sqrt(histmat));
i=i+1;cl(i,:)=clim;
shading interp
xlabel('Range, cm')
ylabel('Speed, cm/s')
title('cricket moving, Dark')
hold on;

nexttile
histmat=histcounts2(mouse_spd(logical(dark.*cricket_moving)), az(logical(dark.*cricket_moving)), Speededges, Azimuthedges);
pcolor(histx(Azimuthedges),Speededges,sqrt(histmat));
i=i+1;cl(i,:)=clim;
shading interp
xlabel('Azimuth, degrees')
ylabel('Speed, cm/s')
title('cricket moving, Dark')
hold on;

nexttile
histmat=histcounts2(az(logical(dark.*cricket_still)), range(logical(dark.*cricket_still)), Azimuthedges, Rangeedges);
pcolor(histx(Azimuthedges),histx(Rangeedges),sqrt(histmat'));
i=i+1;cl(i,:)=clim;
shading interp
xlabel('Azimuth, degrees')
ylabel('Range, cm')
title('cricket still, Dark')
hold on;

nexttile
histmat=histcounts2(mouse_spd(logical(dark.*cricket_still)), range(logical(dark.*cricket_still)), Speededges, Rangeedges);
pcolor(histx(Rangeedges),Speededges,sqrt(histmat));
i=i+1;cl(i,:)=clim;
shading interp
xlabel('Range, cm')
ylabel('Speed, cm/s')
title('cricket still, Dark')
hold on;

nexttile
histmat=histcounts2(mouse_spd(logical(dark.*cricket_still)), az(logical(dark.*cricket_still)), Speededges, Azimuthedges);
pcolor(histx(Azimuthedges),Speededges,sqrt(histmat));
i=i+1;cl(i,:)=clim;
shading interp
xlabel('Azimuth, degrees')
ylabel('Speed, cm/s')
title('cricket still, Dark')
hold on;

nexttile
histmat=histcounts2(az(logical(light.*cricket_moving)), range(logical(light.*cricket_moving)), Azimuthedges, Rangeedges);
pcolor(histx(Azimuthedges),histx(Rangeedges),sqrt(histmat'));
i=i+1;cl(i,:)=clim;
shading interp
xlabel('Azimuth, degrees')
ylabel('Range, cm')
title('cricket moving, Light')
hold on;

nexttile
histmat=histcounts2(mouse_spd(logical(light.*cricket_moving)), range(logical(light.*cricket_moving)), Speededges, Rangeedges);
pcolor(histx(Rangeedges),Speededges,sqrt(histmat));
i=i+1;cl(i,:)=clim;
shading interp
xlabel('Range, cm')
ylabel('Speed, cm/s')
title('cricket moving, Light')
hold on;

nexttile
histmat=histcounts2(mouse_spd(logical(light.*cricket_moving)), az(logical(light.*cricket_moving)), Speededges, Azimuthedges);
pcolor(histx(Azimuthedges),Speededges,sqrt(histmat));
i=i+1;cl(i,:)=clim;
shading interp
xlabel('Azimuth, degrees')
ylabel('Speed, cm/s')
title('cricket moving, Light')
hold on;

nexttile
histmat=histcounts2(az(logical(light.*cricket_still)), range(logical(light.*cricket_still)), Azimuthedges, Rangeedges);
pcolor(histx(Azimuthedges),histx(Rangeedges),sqrt(histmat'));
i=i+1;cl(i,:)=clim;
shading interp
xlabel('Azimuth, degrees')
ylabel('Range, cm')
title('cricket still, Light')
hold on;

nexttile
histmat=histcounts2(mouse_spd(logical(light.*cricket_still)), range(logical(light.*cricket_still)), Speededges, Rangeedges);
pcolor(histx(Rangeedges),histx(Speededges),sqrt(histmat));
i=i+1;cl(i,:)=clim;
shading interp
xlabel('Range, cm')
ylabel('Speed, cm/s')
title('cricket still, Light')
hold on;

nexttile
histmat=histcounts2(mouse_spd(logical(light.*cricket_still)), az(logical(light.*cricket_still)), Speededges, Azimuthedges);
pcolor(histx(Azimuthedges),histx(Speededges),sqrt(histmat));
i=i+1;cl(i,:)=clim;
shading interp
xlabel('Azimuth, degrees')
ylabel('Speed, cm/s')
title('cricket still, Light')
hold on;
shg

for i=1:prod(tl.GridSize),
    nexttile(i)
    set(gca, 'fontsize', 18)
end

if 0 %fix color scale across all plots
    cmax=mean(cl(:,2));
    for i=1:prod(tl.GridSize), nexttile(i), clim([0 cmax]), end
    for i=tl.GridSize(2):tl.GridSize(2):prod(tl.GridSize), nexttile(i), colorbar, end
else
    for i=1:prod(tl.GridSize), nexttile(i), colorbar, end
end


%%%%%%%%%%%%%%%%%%%%%%%%

%%  plot cricket moving contour plots

figure

tl=tiledlayout(2,3, "TileSpacing","compact") ;
set(gcf, "Position", [440 880 1300 840])
title(tl, '(no events or states)')
nexttile

n=2;
lw2=2;
lw=2;

histmat=histcounts2(az(logical(dark.*cricket_moving)), range(logical(dark.*cricket_moving)), Azimuthedges, Rangeedges);
[~,c]=contour(Azimuthedges,Rangeedges,sqrt(histmat'), n, 'm', 'LineWidth', lw);
% L=c.LevelList; %(this snippet is to use a fixed color scale across plots)
% set(c, 'visible', 'off')
% hold on
% contour(Azimuthedges,Rangeedges,sqrt(histmat'), L(1)*[1 1],  'm', 'LineWidth', lw);
% contour(Azimuthedges,Rangeedges,sqrt(histmat'), L(2)*[1 1],  'm', 'LineWidth', lw2);
xlabel('Azimuth, degrees')
ylabel('Range, cm')
hold on
histmat=histcounts2(az(logical(dark.*cricket_still)), range(logical(dark.*cricket_still)), Azimuthedges, Rangeedges);
[~,c]=contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'b', 'LineWidth', lw);
title('dark')
% set(c, 'visible', 'off')
% contour(Azimuthedges,Rangeedges,sqrt(histmat'), c.LevelList(1)*[1 1],  'b', 'LineWidth', lw);
% contour(Azimuthedges,Rangeedges,sqrt(histmat'), c.LevelList(2)*[1 1],  'b', 'LineWidth', lw2);

nexttile
histmat=histcounts2(mouse_spd(logical(dark.*cricket_moving)), range(logical(dark.*cricket_moving)), Speededges, Rangeedges);
contour(Rangeedges,Speededges,sqrt(histmat), n, 'm', 'LineWidth', lw);
xlabel('Range, cm')
ylabel('Speed, cm/s')
hold on
histmat=histcounts2(mouse_spd(logical(dark.*cricket_still)), range(logical(dark.*cricket_still)), Speededges, Rangeedges);
contour(Rangeedges,Speededges,sqrt(histmat), n,  'b', 'LineWidth', lw);
title('dark')

nexttile
histmat=histcounts2(mouse_spd(logical(dark.*cricket_moving)), az(logical(dark.*cricket_moving)), Speededges, Azimuthedges);
contour(Azimuthedges,Speededges,sqrt(histmat), n, 'm', 'LineWidth', lw);
xlabel('Azimuth, degrees')
ylabel('Speed, cm/s')
hold on
histmat=histcounts2(mouse_spd(logical(dark.*cricket_still)), az(logical(dark.*cricket_still)),Speededges,  Azimuthedges);
contour(Azimuthedges,Speededges,sqrt(histmat), n,  'b', 'LineWidth', lw);
legend('cricket moving','cricket still')
title('dark')

nexttile
histmat=histcounts2(az(logical(light.*cricket_moving)), range(logical(light.*cricket_moving)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n, 'm', 'LineWidth', lw);
xlabel('Azimuth, degrees')
ylabel('Range, cm')
hold on
histmat=histcounts2(az(logical(light.*cricket_still)), range(logical(light.*cricket_still)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'b', 'LineWidth', lw);
title('light')

nexttile
histmat=histcounts2(mouse_spd(logical(light.*cricket_moving)), range(logical(light.*cricket_moving)), Speededges, Rangeedges);
contour(Rangeedges,Speededges,sqrt(histmat), n, 'm', 'LineWidth', lw);
xlabel('Range, cm')
ylabel('Speed, cm/s')
hold on
histmat=histcounts2(mouse_spd(logical(light.*cricket_still)), range(logical(light.*cricket_still)), Speededges, Rangeedges);
contour(Rangeedges,Speededges,sqrt(histmat), n,  'b', 'LineWidth', lw);
title('light')

nexttile
histmat=histcounts2(mouse_spd(logical(light.*cricket_moving)), az(logical(light.*cricket_moving)), Speededges, Azimuthedges);
contour(Azimuthedges, Speededges,sqrt(histmat), n, 'm', 'LineWidth', lw);
xlabel('Azimuth, degrees')
ylabel('Speed, cm/s')
hold on
histmat=histcounts2(mouse_spd(logical(light.*cricket_still)), az(logical(light.*cricket_still)),Speededges,  Azimuthedges);
contour(Azimuthedges,Speededges,sqrt(histmat), n,  'b', 'LineWidth', lw);
legend('cricket moving','cricket still')
title('light')

shg
%%%%%%%%%%%%%%%%%%%%%%%%
% more cricket moving contour plots


figure

tl= tiledlayout(2,3, "TileSpacing","compact") ;
set(gcf, "Position", [440 880 1300 840])
title(tl, 'combined dark+light')

nexttile

n=1;
lw=2;

histmat=histcounts2(az(logical(cricket_moving.*chase)), range(logical(cricket_moving.*chase)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n, 'm', 'LineWidth', lw);
xlabel('Azimuth, degrees')
ylabel('Range, cm')
hold on
histmat=histcounts2(az(logical(cricket_still.*chase)), range(logical(cricket_still.*chase)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n, 'b', 'LineWidth', lw);
legend('cricket moving','cricket still')
title('chase')

nexttile
histmat=histcounts2(az(logical(cricket_moving.*approach)), range(logical(cricket_moving.*approach)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'm', 'LineWidth', lw);
xlabel('Azimuth, degrees')
ylabel('Range, cm')
hold on
histmat=histcounts2(az(logical(cricket_still.*approach)), range(logical(cricket_still.*approach)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'b', 'LineWidth', lw);
legend('cricket moving','cricket still')
title('approach')

nexttile
histmat=histcounts2(az(logical(cricket_moving.*pause)), range(logical(cricket_moving.*pause)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'm', 'LineWidth', lw);
xlabel('Azimuth, degrees')
ylabel('Range, cm')
hold on
histmat=histcounts2(az(logical(cricket_still.*pause)), range(logical(cricket_still.*pause)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'b', 'LineWidth', lw);
legend('cricket moving','cricket still')
title('pause')

nexttile
histmat=histcounts2(az(logical(cricket_moving.*wander)), range(logical(cricket_moving.*wander)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'm', 'LineWidth', lw);
xlabel('Azimuth, degrees')
ylabel('Range, cm')
hold on
histmat=histcounts2(az(logical(cricket_still.*wander)), range(logical(cricket_still.*wander)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'b', 'LineWidth', lw);
legend('cricket moving','cricket still')
title('wander')

nexttile
histmat=histcounts2(az(logical(cricket_moving(rangemin_event_frames))), range(logical(cricket_moving(rangemin_event_frames))), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'm', 'LineWidth', lw);
xlabel('Azimuth, degrees')
ylabel('Range, cm')
hold on
histmat=histcounts2(az(logical(cricket_still(rangemin_event_frames))), range(logical(cricket_still(rangemin_event_frames))), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'b', 'LineWidth', lw);
legend('cricket moving','cricket still')
title('rangemin')

nexttile
histmat=histcounts2(az(logical(cricket_moving(cricket_jump_event_frames))), range(logical(cricket_moving(cricket_jump_event_frames))), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'm', 'LineWidth', lw);
xlabel('Azimuth, degrees')
ylabel('Range, cm')
hold on
histmat=histcounts2(az(logical(cricket_still(cricket_jump_event_frames))), range(logical(cricket_still(cricket_jump_event_frames))), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'b', 'LineWidth', lw);
legend('cricket moving','cricket still')
title('cricket jump')

%%%%%%%%%%%%%%%%%%

%% cricket moving 2D scatterplot of random subsets of state frames
% randomized plotting order (with color matrix) so plotting order doesn't obscure

numpoints=10000;
alph=.15;


x1=az(logical(dark.*cricket_moving));
y1=range(logical(dark.*cricket_moving));
z1=mouse_spd(logical(dark.*cricket_moving));

x2=az(logical(dark.*cricket_still));
y2=range(logical(dark.*cricket_still));
z2=mouse_spd(logical(dark.*cricket_still));


%random subset
random_idx = randperm(length(x1), numpoints);
x1 = x1(random_idx);
y1 = y1(random_idx);
z1 = z1(random_idx);
random_idx = randperm(length(x2), numpoints);
x2 = x2(random_idx);
y2 = y2(random_idx);
z2 = z2(random_idx);



%combine and randomize plotting order
X=[x1; x2]; % az
Y=[y1; y2]; % range
Z=[z1; z2]; % speed
shuffle_idx = randperm(length(X));

% Create a Matching Color Matrix (RGB Triplets) ---
% Red is [1 0 0], Blue is [0 0 1]
colors_red  = repmat([1, 0, 0], numpoints, 1);
colors_cyan = repmat([0, 1, 1], numpoints, 1);
colors_green = repmat([0, 1, 0], numpoints, 1);
colors_mag = repmat([1, 0, 1], numpoints, 1);
colors_blue = repmat([0, 0, 1], numpoints, 1);
colors_orange = repmat([1, .65, 0], numpoints, 1);

C = [colors_mag; colors_blue];

figure
tiledlayout(1,3, "TileSpacing","compact")
set(gcf, "Position", [440 880 1300 420])

nexttile
h=scatter(X(shuffle_idx), Y(shuffle_idx), 20, C(shuffle_idx, :), 'filled', 'MarkerFaceAlpha', alph);
title(sprintf('dark, %d random frames', numpoints))
xlabel('Azimuth, degrees')
ylabel('Range, cm')
xlim([0 180]) %az
ylim([0 60]) %range
set(gca, 'fontsize', 18)

nexttile
h=scatter(Y(shuffle_idx), Z(shuffle_idx), 20, C(shuffle_idx, :), 'filled', 'MarkerFaceAlpha', alph);
title(sprintf('dark, %d random frames', numpoints))
xlabel('Range, cm')
ylabel('Speed, cm/s')
xlim([0 60]) %range
ylim([0 60]) %speed
set(gca, 'fontsize', 18)

nexttile
h=scatter(X(shuffle_idx), Z(shuffle_idx), 20, C(shuffle_idx, :), 'filled', 'MarkerFaceAlpha', alph);
title(sprintf('dark, %d random frames', numpoints))
xlabel('Azimuth, degrees')
ylabel('Speed, cm/s')
xlim([0 180]) %az
ylim([0 60]) %speed
set(gca, 'fontsize', 18)

% 1. Create invisible dummy plots with the exact same visual properties
hold on
h_mag = scatter(NaN, NaN, 20, [1, 0, 1], 'filled');
h_blue = scatter(NaN, NaN, 20, [0, 0, 1], 'filled');

% 2. Pass only these fake handles to the legend function
legend([h_mag, h_blue], {'cricket moving', 'cricket still'}, 'Location', 'best');

xlim([0 60]) %speed



%% cricket moving 3D scatterplot of random subsets of state frames

cd '/Users/wehr/Documents/Projects/Molly Shallow'
figure
h=scatter3(X(shuffle_idx), Y(shuffle_idx), Z(shuffle_idx), 20, C(shuffle_idx, :), 'filled', 'MarkerFaceAlpha', alph);
title(sprintf('dark, %d random frames', numpoints))
xlabel('Azimuth, degrees')
ylabel('Range, cm')
zlabel('Speed, cm/s')


% 1. Create invisible dummy plots with the exact same visual properties
hold on
h_mag = scatter3(NaN, NaN, NaN, 20, [1, 0, 1], 'filled');
h_blue = scatter3(NaN, NaN, NaN, 20, [0, 0, 1], 'filled');

% 2. Pass only these fake handles to the legend function
legend([h_mag, h_blue], {'cricket moving', 'cricket still'}, 'Location', 'best');

xlim([0 180])%az
ylim([0 60])%range
zlim([0 60])%spd
xticks([0:30:180])
yticks([0:20:60])
zticks([0:20:60])

axis vis3d;  % <-- MAGIC LINE: Freezes the aspect ratio for 3D rotation

% --- 2. Initialize the Video Writer ---
% 'MPEG-4' creates an MP4 file which works everywhere (Windows, Mac, Mobile)
if 0 %create movie
    video_filename = '3d_cricketmoving.mp4';
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 30;  % 30 frames per second for smooth playback
    v.Quality = 95;    % High quality compression (out of 100)
    open(v);

    elevation_angle = 30;
    azimuth_angles = 45:1:720+45; % Adjust the step (middle number) to change speed
    % --- 4. Loop, Update View, and Capture Frames ---
    for azi = azimuth_angles
        if azi>360
            elevation_angle = 30 + 30 * sin(deg2rad(azi));so I
        end

        % Update the camera view angle
        view(azi, elevation_angle);
        fprintf('\n%d', azi')
        % Force MATLAB to draw the update immediately
        drawnow;

        % Capture the entire figure window as a frame
        frame = getframe(gcf);

        % Write the frame to the video file
        writeVideo(v, frame);
    end
    close(v);
    disp(['Movie saved successfully as: ', video_filename]);
end
if 1 %create gif version of movie
 
    gif_filename = '3d_cricketmoving.gif';
    delay_time = 1/30; % Matches a 30 FPS playback speed

    elevation_angle = 30;
    azimuth_angles = 45:1:720+45; % Adjust the step (middle number) to change speed
    for i = 1:length(azimuth_angles)
        if azimuth_angles(i)>360
            elevation_angle = 30 + 30 * sin(deg2rad(azimuth_angles(i)));
        end
        view(azimuth_angles(i), elevation_angle);
        drawnow;

        % Capture and convert the frame to an image matrix
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256); % Convert to indexed color map
        fprintf('\n%d', i')

        % Write to the GIF file
        if i == 1
            % 'Loopcount', inf forces it to loop forever natively everywhere
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
        end
    end
    % you can reduce gif filesize (e.g. to send with discord or put in google slides) like this:
    %ffmpeg -i 3d_cricketmoving.gif -vf "fps=15,scale=600:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" 3d_cricketmoving_small.gif
end


%%%%%%%%%%
%% how many state frames are there in cricket-moving vs cricket still?


figure
tiledlayout(1,2, "TileSpacing","compact")
set(gcf, "Position", [440 880 1300 570])
nexttile

hold on
h1=scatter(sum(cricket_moving.*chase.*dark), sum(cricket_still.*chase.*dark), 'r', 'filled');
h2=scatter(sum(cricket_moving.*approach.*dark), sum(cricket_still.*approach.*dark), 'm', 'filled');
h3=scatter(sum(cricket_moving(dark(rangemin_event_frames))), sum(cricket_still(dark(rangemin_event_frames))), 'b', 'filled');
h4=scatter(sum(cricket_moving.*wander.*dark), sum(cricket_still.*wander.*dark), 'g', 'filled');
h5=scatter(sum(cricket_moving.*pause.*dark), sum(cricket_still.*pause.*dark), 'c', 'filled');
h6=scatter(sum(cricket_moving(dark(cricket_jump_event_frames))), sum(cricket_still(dark(cricket_jump_event_frames))), 'y','filled');
set(h6, 'markerfacecolor', [1, .65, 0])
set([h1 h2 h3 h4 h5 h6], 'SizeData', 200)
xlabel('cricket moving')
ylabel('cricket still')
xl=xlim;
yl=ylim;
maxl=max([xl, yl]);
xlim([0 maxl])
ylim([0 maxl])
line([0 maxl], [0 maxl])
% legend([h1 h2 h3 h4 h5 h6], {'chase', 'approach', 'rangemin', 'wander', 'pause', 'cricketjump'}, 'Location', 'best');
title('number of frames, dark')
set(gca, 'fontsize', 16)

nexttile
hold on
h1=scatter(sum(cricket_moving.*chase.*light), sum(cricket_still.*chase.*light), 'r', 'filled');
h2=scatter(sum(cricket_moving.*approach.*light), sum(cricket_still.*approach.*light), 'm', 'filled');
h3=scatter(sum(cricket_moving(light(rangemin_event_frames))), sum(cricket_still(light(rangemin_event_frames))), 'b', 'filled');
h4=scatter(sum(cricket_moving.*wander.*light), sum(cricket_still.*wander.*light), 'g', 'filled');
h5=scatter(sum(cricket_moving.*pause.*light), sum(cricket_still.*pause.*light), 'c', 'filled');
h6=scatter(sum(cricket_moving(light(cricket_jump_event_frames))), sum(cricket_still(light(cricket_jump_event_frames))), 'y','filled');
set(h6, 'markerfacecolor', [1, .65, 0])
set([h1 h2 h3 h4 h5 h6], 'SizeData', 200)
xlabel('cricket moving')
ylabel('cricket still')
xl=xlim;
yl=ylim;
maxl=max([xl, yl]);
xlim([0 maxl])
ylim([0 maxl])
line([0 maxl], [0 maxl])
legend([h1 h2 h3 h4 h5 h6], {'chase', 'approach', 'rangemin', 'wander', 'pause', 'cricketjump'}, 'Location', 'best');
title('number of frames, light')
set(gca, 'fontsize', 16)

%zoom in on cricketjump and rangemin
nexttile(1)
xlim([0 5e3])
ylim([0 5e3])
set(gca, 'fontsize', 32)
nexttile(2)
xlim([0 5e3])
ylim([0 5e3])
set(gca, 'fontsize', 32)

%bar graph of how many frames are in cricket moving vs still
categories_list = {'cricket moving dark', 'cricket still dark', 'cricket moving light', 'cricket still light'};
values = [sum(cricket_moving.*dark), sum(cricket_still.*dark), sum(cricket_moving.*light), sum(cricket_still.*light)];
% Convert to categorical
desired_order = categories_list;
cat_data = categorical(categories_list, desired_order);
figure;
b = bar(cat_data, values, 'FaceColor', 'flat');
ylabel('frame count');
my_colors = [0.55, 0.00, 0.00;    % dark red
    0.00, 0.00, 0.55;    % dark blue
    1 0 0;    % bright red
    0 1 1]; % cyan
% Assign the colors to the bar object
b.CData = my_colors;


%%%%%%%%%%%%
%%  2D joint probability heat maps
% %loop through all states

%statenames={'chase', 'approach', 'wander', 'pause'};

thresholds.chase.az=[0 30];
thresholds.chase.mouse_spd=[10 60];
thresholds.chase.range=[5 10];

thresholds.approach.az=[0 30];
thresholds.approach.mouse_spd=[5 60];
thresholds.approach.range=[0 60];  %wide open

thresholds.wander.az=[30 180];
thresholds.wander.mouse_spd=[5 15];
thresholds.wander.range=[15 60];

thresholds.pause.az=[0 180]; %wide open
thresholds.pause.mouse_spd=[0 5];
thresholds.pause.range=[0 60]; %wide open

thresholds.stalk.az=[30 180]; %wide open
thresholds.stalk.mouse_spd=[3 15];
thresholds.stalk.range=[5 40]; %wide open

thresholds.none.az=[0 0];
thresholds.none.range=[0 0];
thresholds.none.mouse_spd=[0 0];

% dummy thresholds
thresholds.hotpursuit.az=[0 0];
thresholds.hotpursuit.range=[0 0];
thresholds.hotpursuit.mouse_spd=[0 0];
thresholds.follow=thresholds.hotpursuit;

for statename=statenames

    state=eval(statename{:});
    dosqrt=1;

    figure
    tiledlayout(4,3, "TileSpacing","compact")
    set(gcf, "Position", [440 880 1300 4*420])



    %dark & state
    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(az(dark&state), range(dark&state), Azimuthedges, Rangeedges));
    else
        histmat=(histcounts2(az(dark&state), range(dark&state), Azimuthedges, Rangeedges));
    end
    pcolor(histx(Azimuthedges),histx(Rangeedges),histmat');
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Range, cm')
    title(sprintf('%s dark', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'az');
    y=getfield(thresholds, statename{:}, 'range');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');

    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(mouse_spd(dark&state), range(dark&state), Speededges, Rangeedges));
    else
        histmat=(histcounts2(mouse_spd(dark&state), range(dark&state), Speededges, Rangeedges));
    end
    pcolor(histx(Rangeedges),histx(Speededges), histmat);
    shading interp
    xlabel('Range, cm')
    ylabel('Speed, cm/s')
    title(sprintf('%s dark', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'range');
    y=getfield(thresholds, statename{:}, 'mouse_spd');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');

    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(mouse_spd(dark&state), az(dark&state), Speededges, Azimuthedges));
    else
        histmat=(histcounts2(mouse_spd(dark&state), az(dark&state), Speededges, Azimuthedges));
    end
    pcolor(histx(Azimuthedges),histx(Speededges),histmat);
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Speed, cm/s')
    title(sprintf('%s dark', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'az');
    y=getfield(thresholds, statename{:}, 'mouse_spd');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');
    h=colorbar;
    ylabel(h, sprintf('sqrt=%d', dosqrt));

    %dark & ~state
    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(az(dark&~state), range(dark&~state), Azimuthedges, Rangeedges));
    else
        histmat=(histcounts2(az(dark&~state), range(dark&~state), Azimuthedges, Rangeedges));
    end
    pcolor(histx(Azimuthedges),histx(Rangeedges),histmat');
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Range, cm')
    title(sprintf('%s FALSE, dark', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'az');
    y=getfield(thresholds, statename{:}, 'range');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');

    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(mouse_spd(dark&~state), range(dark&~state), Speededges, Rangeedges));
    else
        histmat=(histcounts2(mouse_spd(dark&~state), range(dark&~state), Speededges, Rangeedges));
    end
    pcolor(histx(Rangeedges),histx(Speededges), histmat);
    shading interp
    xlabel('Range, cm')
    ylabel('Speed, cm/s')
    title(sprintf('%s FALSE, dark', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'range');
    y=getfield(thresholds, statename{:}, 'mouse_spd');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');

    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(mouse_spd(dark&~state), az(dark&~state), Speededges, Azimuthedges));
    else
        histmat=(histcounts2(mouse_spd(dark&~state), az(dark&~state), Speededges, Azimuthedges));
    end
    pcolor(histx(Azimuthedges),histx(Speededges),histmat);
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Speed, cm/s')
    title(sprintf('%s FALSE, dark', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'az');
    y=getfield(thresholds, statename{:}, 'mouse_spd');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');
    h=colorbar;
    ylabel(h, sprintf('sqrt=%d', dosqrt));

    %light & state
    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(az(light&state), range(light&state), Azimuthedges, Rangeedges));
    else
        histmat=(histcounts2(az(light&state), range(light&state), Azimuthedges, Rangeedges));
    end
    pcolor(histx(Azimuthedges),histx(Rangeedges),histmat');
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Range, cm')
    title(sprintf('%s light', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'az');
    y=getfield(thresholds, statename{:}, 'range');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');

    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(mouse_spd(light&state), range(light&state), Speededges, Rangeedges));
    else
        histmat=(histcounts2(mouse_spd(light&state), range(light&state), Speededges, Rangeedges));
    end
    pcolor(histx(Rangeedges),histx(Speededges), histmat);
    shading interp
    xlabel('Range, cm')
    ylabel('Speed, cm/s')
    title(sprintf('%s light', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'range');
    y=getfield(thresholds, statename{:}, 'mouse_spd');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');

    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(mouse_spd(light&state), az(light&state), Speededges, Azimuthedges));
    else
        histmat=(histcounts2(mouse_spd(light&state), az(light&state), Speededges, Azimuthedges));
    end
    pcolor(histx(Azimuthedges),histx(Speededges),histmat);
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Speed, cm/s')
    title(sprintf('%s light', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'az');
    y=getfield(thresholds, statename{:}, 'mouse_spd');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');
    h=colorbar;
    ylabel(h, sprintf('sqrt=%d', dosqrt));

    %light & ~state
    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(az(light&~state), range(light&~state), Azimuthedges, Rangeedges));
    else
        histmat=(histcounts2(az(light&~state), range(light&~state), Azimuthedges, Rangeedges));
    end
    pcolor(histx(Azimuthedges),histx(Rangeedges),histmat');
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Range, cm')
    title(sprintf('%s FALSE, light', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'az');
    y=getfield(thresholds, statename{:}, 'range');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');

    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(mouse_spd(light&~state), range(light&~state), Speededges, Rangeedges));
    else
        histmat=(histcounts2(mouse_spd(light&~state), range(light&~state), Speededges, Rangeedges));
    end
    pcolor(histx(Rangeedges),histx(Speededges), histmat);
    shading interp
    xlabel('Range, cm')
    ylabel('Speed, cm/s')
    title(sprintf('%s FALSE, light', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'range');
    y=getfield(thresholds, statename{:}, 'mouse_spd');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');

    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(mouse_spd(light&~state), az(light&~state), Speededges, Azimuthedges));
    else
        histmat=(histcounts2(mouse_spd(light&~state), az(light&~state), Speededges, Azimuthedges));
    end
    pcolor(histx(Azimuthedges),histx(Speededges),histmat);
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Speed, cm/s')
    title(sprintf('%s FALSE, light', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'az');
    y=getfield(thresholds, statename{:}, 'mouse_spd');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');
    h=colorbar;
    ylabel(h, sprintf('sqrt=%d', dosqrt));

end
set(get((get(gcf, 'Children')), 'Children'), 'fontsize', 18)

%%%%%%%%%%%%
%%  more 2D joint probability heat maps
% %loop through all events


statenames={'rangemin_event_frames', 'cricket_jump_event_frames'}
for statename=statenames

    state=eval(statename{:});
    dosqrt=0;

    figure
    tiledlayout(2,3, "TileSpacing","compact")
    set(gcf, "Position", [440 880 1300 2*420])



    %dark & state
    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(az(dark(state)), range(dark(state)), Azimuthedges, Rangeedges));
    else
        histmat=(histcounts2(az(dark(state)), range(dark(state)), Azimuthedges, Rangeedges));
    end
    pcolor(histx(Azimuthedges),histx(Rangeedges),histmat');
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Range, cm')
    title(sprintf('%s dark', erase(statename{:}, '_event_frames')), 'interpreter', 'none')

    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(mouse_spd(dark(state)), range(dark(state)), Speededges, Rangeedges));
    else
        histmat=(histcounts2(mouse_spd(dark(state)), range(dark(state)), Speededges, Rangeedges));
    end
    pcolor(histx(Rangeedges),histx(Speededges),histmat);
    shading interp
    xlabel('Range, cm')
    ylabel('Speed, cm/s')
    title(sprintf('%s dark', erase(statename{:}, '_event_frames')), 'interpreter', 'none')

    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(mouse_spd(dark(state)), az(dark(state)), Speededges, Azimuthedges));
    else
        histmat=(histcounts2(mouse_spd(dark(state)), az(dark(state)), Speededges, Azimuthedges));
    end
    pcolor(histx(Azimuthedges),histx(Speededges),histmat);
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Speed, cm/s')
    title(sprintf('%s dark', erase(statename{:}, '_event_frames')), 'interpreter', 'none')
    h=colorbar;
    ylabel(h, sprintf('sqrt=%d', dosqrt));


    %light & state
    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(az(light(state)), range(light(state)), Azimuthedges, Rangeedges));
    else
        histmat=(histcounts2(az(light(state)), range(light(state)), Azimuthedges, Rangeedges));
    end
    pcolor(histx(Azimuthedges),histx(Rangeedges),histmat');
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Range, cm')
    title(sprintf('%s light', erase(statename{:}, '_event_frames')), 'interpreter', 'none')

    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(mouse_spd(light(state)), range(light(state)), Speededges, Rangeedges));
    else
        histmat=(histcounts2(mouse_spd(light(state)), range(light(state)), Speededges, Rangeedges));
    end
    pcolor(histx(Rangeedges),histx(Speededges),histmat);
    shading interp
    xlabel('Range, cm')
    ylabel('Speed, cm/s')
    title(sprintf('%s light', erase(statename{:}, '_event_frames')), 'interpreter', 'none')

    nexttile
    if dosqrt
        histmat=sqrt(histcounts2(mouse_spd(light(state)), az(light(state)), Speededges, Azimuthedges));
    else
        histmat=(histcounts2(mouse_spd(light(state)), az(light(state)), Speededges, Azimuthedges));
    end
    pcolor(histx(Azimuthedges),histx(Speededges),histmat);
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Speed, cm/s')
    title(sprintf('%s light', erase(statename{:}, '_event_frames')), 'interpreter', 'none')
    h=colorbar;
    ylabel(h, sprintf('sqrt=%d', dosqrt));

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% contour plots
cd ~

figure

tiledlayout(2,3, "TileSpacing","compact")
set(gcf, "Position", [440 880 1300 840])
nexttile

n=1;
lw=2;

histmat=histcounts2(az(logical(dark.*chase)), range(logical(dark.*chase)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n, 'r', 'LineWidth', lw);
xlabel('Azimuth, degrees')
ylabel('Range, cm')
hold on
histmat=histcounts2(az(logical(dark.*approach)), range(logical(dark.*approach)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'm', 'LineWidth', lw);
histmat=histcounts2(az(logical(dark.*pause)), range(logical(dark.*pause)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'b', 'LineWidth', lw);
histmat=histcounts2(az(logical(dark.*wander)), range(logical(dark.*wander)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'g', 'LineWidth', lw);
histmat=histcounts2(az(logical(dark(rangemin_event_frames))), range(logical(dark(rangemin_event_frames))), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'c', 'LineWidth', lw);
legend('chase', 'approach', 'pause', 'wander',  'rangemin')
title('dark')

nexttile
histmat=histcounts2(mouse_spd(logical(dark.*chase)), range(logical(dark.*chase)), Speededges, Rangeedges);
contour(Rangeedges,Speededges, sqrt(histmat), n, 'r', 'LineWidth', lw);
xlabel('Range, cm')
ylabel('Speed, cm/s')
hold on
histmat=histcounts2(mouse_spd(logical(dark.*approach)), range(logical(dark.*approach)), Speededges, Rangeedges);
contour(Rangeedges,Speededges,sqrt(histmat) , n,  'm', 'LineWidth', lw);
histmat=histcounts2(mouse_spd(logical(dark.*pause)), range(logical(dark.*pause)), Speededges, Rangeedges);
contour(Rangeedges,Speededges,sqrt(histmat) , n,  'b', 'LineWidth', lw);
histmat=histcounts2(mouse_spd(logical(dark.*wander)), range(logical(dark.*wander)), Speededges, Rangeedges);
contour(Rangeedges, Speededges, sqrt(histmat), n,  'g', 'LineWidth', lw);
histmat=histcounts2(mouse_spd(logical(dark(rangemin_event_frames))), range(logical(dark(rangemin_event_frames))), Speededges, Rangeedges);
contour(Rangeedges,Speededges,sqrt(histmat) , n,  'c', 'LineWidth', lw);
legend('chase', 'approach', 'pause', 'wander',  'rangemin')
title('dark')

nexttile
histmat=histcounts2(mouse_spd(logical(dark.*chase)), az(logical(dark.*chase)), Speededges, Azimuthedges);
contour(Azimuthedges,Speededges,sqrt(histmat), n, 'r', 'LineWidth', lw);
xlabel('Azimuth, degrees')
ylabel('Speed, cm/s')
hold on
histmat=histcounts2(mouse_spd(logical(dark.*approach)), az(logical(dark.*approach)),Speededges,  Azimuthedges);
contour(Rangeedges,Speededges,sqrt(histmat) , n,  'm', 'LineWidth', lw);
histmat=histcounts2(mouse_spd(logical(dark.*pause)), az(logical(dark.*pause)), Speededges, Azimuthedges);
contour(Azimuthedges,Speededges,sqrt(histmat), n,  'b', 'LineWidth', lw);
histmat=histcounts2(mouse_spd(logical(dark.*wander)), az(logical(dark.*wander)), Speededges, Azimuthedges);
contour(Azimuthedges,Speededges,sqrt(histmat), n,  'g', 'LineWidth', lw);
histmat=histcounts2(mouse_spd(logical(dark(rangemin_event_frames))), az(logical(dark(rangemin_event_frames))), Speededges, Azimuthedges);
contour(Rangeedges,Speededges,sqrt(histmat) , n,  'c', 'LineWidth', lw);
legend('chase', 'approach', 'pause', 'wander',  'rangemin')
title('dark')

nexttile
histmat=histcounts2(az(logical(light.*chase)), range(logical(light.*chase)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n, 'r', 'LineWidth', lw);
xlabel('Azimuth, degrees')
ylabel('Range, cm')
hold on
histmat=histcounts2(az(logical(light.*approach)), range(logical(light.*approach)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'm', 'LineWidth', lw);
histmat=histcounts2(az(logical(light.*pause)), range(logical(light.*pause)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'b', 'LineWidth', lw);
histmat=histcounts2(az(logical(light.*wander)), range(logical(light.*wander)), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'g', 'LineWidth', lw);
histmat=histcounts2(az(logical(light(rangemin_event_frames))), range(logical(light(rangemin_event_frames))), Azimuthedges, Rangeedges);
contour(Azimuthedges,Rangeedges,sqrt(histmat'), n,  'c', 'LineWidth', lw);legend('chase', 'approach', 'pause', 'wander',  'rangemin')
title('light')

nexttile
histmat=histcounts2(mouse_spd(logical(light.*chase)), range(logical(light.*chase)), Speededges, Rangeedges);
contour(Rangeedges,Speededges,sqrt(histmat) , n, 'r', 'LineWidth', lw);
xlabel('Range, cm')
ylabel('Speed, cm/s')
hold on
histmat=histcounts2(mouse_spd(logical(light.*approach)), range(logical(light.*approach)), Speededges, Rangeedges);
contour(Rangeedges,Speededges,sqrt(histmat) , n,  'm', 'LineWidth', lw);
histmat=histcounts2(mouse_spd(logical(light.*pause)), range(logical(light.*pause)), Speededges, Rangeedges);
contour(Rangeedges,Speededges,sqrt(histmat) , n,  'b', 'LineWidth', lw);
histmat=histcounts2(mouse_spd(logical(light.*wander)), range(logical(light.*wander)), Speededges, Rangeedges);
contour(Rangeedges,Speededges, sqrt(histmat), n,  'g', 'LineWidth', lw);
histmat=histcounts2(mouse_spd(logical(light(rangemin_event_frames))), range(logical(light(rangemin_event_frames))), Speededges, Rangeedges);
contour(Rangeedges,Speededges,sqrt(histmat) , n,  'c', 'LineWidth', lw);
legend('chase', 'approach', 'pause', 'wander',  'rangemin')
title('light')

nexttile
histmat=histcounts2(mouse_spd(logical(light.*chase)), az(logical(light.*chase)), Speededges, Azimuthedges);
contour(Azimuthedges,Speededges,sqrt(histmat), n, 'r', 'LineWidth', lw);
xlabel('Azimuth, degrees')
ylabel('Speed, cm/s')
hold on
histmat=histcounts2(mouse_spd(logical(light.*approach)), az(logical(light.*approach)), Speededges, Azimuthedges);
contour(Azimuthedges,Speededges,sqrt(histmat), n,  'm', 'LineWidth', lw);
histmat=histcounts2(mouse_spd(logical(light.*pause)), az(logical(light.*pause)), Speededges, Azimuthedges);
contour(Azimuthedges,Speededges,sqrt(histmat), n,  'b', 'LineWidth', lw);
histmat=histcounts2(mouse_spd(logical(light.*wander)), az(logical(light.*wander)), Speededges, Azimuthedges);
contour(Azimuthedges,Speededges,sqrt(histmat), n,  'g', 'LineWidth', lw);
histmat=histcounts2(mouse_spd(logical(light(rangemin_event_frames))), az(logical(light(rangemin_event_frames))), Speededges, Azimuthedges);
contour(Rangeedges,Speededges,sqrt(histmat) , n,  'c', 'LineWidth', lw);
legend('chase', 'approach', 'pause', 'wander',  'rangemin')
title('light')

shg

%% 2D scatterplot of random subsets of state frames
% randomized plotting order (with color matrix) so plotting order doesn't obscure
numpoints=10000;
alph=.25;

x1=az(chase);
y1=range(chase);
z1=mouse_spd(chase);

x2=az(pause);
y2=range(pause);
z2=mouse_spd(pause);

x3=az(wander);
y3=range(wander);
z3=mouse_spd(wander);

x4=az(rangemin_event_frames);
y4=range(rangemin_event_frames);
z4=mouse_spd(rangemin_event_frames);

x5=az(approach);
y5=range(approach);
z5=mouse_spd(approach);

x6=az(cricket_jump_event_frames);
y6=range(cricket_jump_event_frames);
z6=mouse_spd(cricket_jump_event_frames);



%random subset
random_idx = randperm(length(x1), numpoints);
x1 = x1(random_idx);
y1 = y1(random_idx);
z1 = z1(random_idx);
random_idx = randperm(length(x2), numpoints);
x2 = x2(random_idx);
y2 = y2(random_idx);
z2 = z2(random_idx);
random_idx = randperm(length(x3), numpoints);
x3 = x3(random_idx);
y3 = y3(random_idx);
z3 = z3(random_idx);
%no x4 because rangemin is not even numpoints long (because it's events not a state)
random_idx = randperm(length(x5), numpoints);
x5 = x5(random_idx);
y5 = y5(random_idx);
z5 = z5(random_idx);
%no x6 because cricketjump is not even numpoints long (because it's events not a state)
% random_idx = randperm(length(x6), numpoints);
% x6 = x6(random_idx);
% y6 = y6(random_idx);
% z6 = z6(random_idx);


%combine and randomize plotting order
X=[x1; x2; x3; x4; x5; x6]; % az
Y=[y1; y2; y3; y4; y5; y6]; % range
Z=[z1; z2; z3; z4; z5; z6]; % speed
shuffle_idx = randperm(length(X));

% Create a Matching Color Matrix (RGB Triplets) ---
% Red is [1 0 0], Blue is [0 0 1]
colors_red  = repmat([1, 0, 0], numpoints, 1);
colors_cyan = repmat([0, 1, 1], numpoints, 1);
colors_green = repmat([0, .8, 0], numpoints, 1);
colors_mag = repmat([1, 0, 1], numpoints, 1);
colors_blue = repmat([0, 0, 1], length(x4), 1);
colors_orange = repmat([1, .65, 0], length(x6), 1);

C = [colors_red; colors_cyan ; colors_green; colors_blue; colors_mag; colors_orange];

figure
tiledlayout(1,3, "TileSpacing","compact")
set(gcf, "Position", [440 880 1300 420])

nexttile
h=scatter(X(shuffle_idx), Y(shuffle_idx), 20, C(shuffle_idx, :), 'filled', 'MarkerFaceAlpha', alph);
title(sprintf('%d random frames', numpoints))
xlabel('Azimuth, degrees')
ylabel('Range, cm')
xlim([0 180]) %az
ylim([0 60]) %range

nexttile
h=scatter(Y(shuffle_idx), Z(shuffle_idx), 20, C(shuffle_idx, :), 'filled', 'MarkerFaceAlpha', alph);
title(sprintf('%d random frames', numpoints))
xlabel('Range, cm')
ylabel('Speed, cm/s')
xlim([0 60]) %speed
ylim([0 60]) %range

nexttile
h=scatter(X(shuffle_idx), Z(shuffle_idx), 20, C(shuffle_idx, :), 'filled', 'MarkerFaceAlpha', alph);
title(sprintf('%d random frames', numpoints))
xlabel('Azimuth, degrees')
ylabel('Speed, cm/s')
xlim([0 60]) %speed
ylim([0 180]) %az

% 1. Create invisible dummy plots with the exact same visual properties
hold on
h_red  = scatter(NaN, NaN, 20, [1, 0, 0], 'filled');
h_blue = scatter(NaN, NaN, 20, [0, 0, 1], 'filled');
h_green = scatter(NaN, NaN, 20, [0, .8, 0], 'filled');
h_cyan = scatter(NaN, NaN, 20, [0, 1, 1], 'filled');
h_mag = scatter(NaN, NaN, 20, [1, 0, 1], 'filled');
h_orange = scatter(NaN, NaN, 20, [1, .65, 0], 'filled');

% 2. Pass only these fake handles to the legend function
legend([h_red, h_mag, h_blue, h_green, h_cyan, h_orange], {'chase', 'approach', 'rangemin', 'wander', 'pause', 'cricketjump'}, 'Location', 'best');

xlim([0 60]) %speed


%%  3D scatterplot of random subsets of state frames

cd '/Users/wehr/Documents/Projects/Molly Shallow'
figure
h=scatter3(X(shuffle_idx), Y(shuffle_idx), Z(shuffle_idx), 20, C(shuffle_idx, :), 'filled', 'MarkerFaceAlpha', alph);
title(sprintf('%d random frames', numpoints))
xlabel('Azimuth, degrees')
ylabel('Range, cm')
zlabel('Speed, cm/s')


% 1. Create invisible dummy plots with the exact same visual properties
hold on
h_red  = scatter3(NaN,NaN, NaN, 20, [1, 0, 0], 'filled');
h_blue = scatter3(NaN, NaN, NaN, 20, [0, 0, 1], 'filled');
h_green = scatter3(NaN, NaN, NaN,20, [0, .8, 0], 'filled');
h_cyan = scatter3(NaN, NaN, NaN, 20, [0, 1, 1], 'filled');
h_mag = scatter3(NaN, NaN, NaN, 20, [1, 0, 1], 'filled');
h_orange = scatter3(NaN, NaN, NaN, 20, [1, .65, 0], 'filled');

% 2. Pass only these fake handles to the legend function
legend([h_red, h_mag, h_blue, h_green, h_cyan, h_orange], {'chase', 'approach', 'rangemin', 'wander', 'pause', 'cricketjump'}, 'Location', 'best');

xlim([0 180])%az
ylim([0 60])%range
zlim([0 60])%spd
xticks([0:30:180])
yticks([0:20:60])
zticks([0:20:60])

axis vis3d;  % <-- MAGIC LINE: Freezes the aspect ratio for 3D rotation

% --- 2. Initialize the Video Writer ---
% 'MPEG-4' creates an MP4 file which works everywhere (Windows, Mac, Mobile)
if 0 %create movie

    video_filename = '3d_spin.mp4';
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 30;  % 30 frames per second for smooth playback
    v.Quality = 95;    % High quality compression (out of 100)
    open(v);

    elevation_angle = 30;
    azimuth_angles = 45:1:720+45; % Adjust the step (middle number) to change speed
    % --- 4. Loop, Update View, and Capture Frames ---
    for azi = azimuth_angles
        if azi>360
            elevation_angle = 30 + 30 * sin(deg2rad(azi));
        end

        % Update the camera view angle
        view(azi, elevation_angle);
        fprintf('\n%d', azi')
        % Force MATLAB to draw the update immediately
        drawnow;

        % Capture the entire figure window as a frame
        frame = getframe(gcf);

        % Write the frame to the video file
        writeVideo(v, frame);
    end
    close(v);
    disp(['Movie saved successfully as: ', video_filename]);
end

% ____________
%% gif allows looping by default, but quicktime won't play, use a browser

if 0 %create movie
    gif_filename = '3d_spin.gif';
    delay_time = 1/30; % Matches a 30 FPS playback speed

    for i = 1:length(azimuth_angles)
        if azimuth_angles(i)>360
            elevation_angle = 30 + 30 * sin(deg2rad(azimuth_angles(i)));
        end
        view(azimuth_angles(i), elevation_angle);
        drawnow;

        % Capture and convert the frame to an image matrix
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256); % Convert to indexed color map
        fprintf('\n%d', i')

        % Write to the GIF file
        if i == 1
            % 'Loopcount', inf forces it to loop forever natively everywhere
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
        end
    end
    % you can reduce gif filesize (e.g. to send with discord) like this:
    %ffmpeg -i 3D-spin.gif -vf "fps=15,scale=600:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" 3D-spin_small.gif
end







%%%%%%%%%%%%
%%  laser on-off
%loop through all states

dosqrt=0;

statenames={'chase', 'approach', 'wander', 'pause'};
for statename=statenames

    state=eval(statename{:});

    figure
    tiledlayout(1,3, "TileSpacing","compact")
    set(gcf, "Position", [440 880 1300 420])
    nexttile
    if dosqrt
        histmaton=sqrt(histcounts2(az(dark&state&laseron), range(dark&state&laseron), Azimuthedges, Rangeedges));
        histmatoff=sqrt(histcounts2(az(dark&state&~laseron), range(dark&state&~laseron), Azimuthedges, Rangeedges));
    else
        histmaton=(histcounts2(az(dark&state&laseron), range(dark&state&laseron), Azimuthedges, Rangeedges));
        histmatoff=(histcounts2(az(dark&state&~laseron), range(dark&state&~laseron), Azimuthedges, Rangeedges));
    end
    pcolor(histx(Azimuthedges),histx(Rangeedges),(histmaton'-histmatoff'));
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Range, cm')
    title(sprintf('%s laser on-off (dark)', statename{:}))
    cl=clim;
    clim([-1 1]*max(abs(cl)))

    hold on;
    x=getfield(thresholds, statename{:}, 'az');
    y=getfield(thresholds, statename{:}, 'range');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');

    nexttile
    if dosqrt
        histmaton=sqrt(histcounts2(mouse_spd(dark&state&laseron), range(dark&state&laseron), Speededges, Rangeedges));
        histmatoff=sqrt(histcounts2(mouse_spd(dark&state&~laseron), range(dark&state&~laseron), Speededges, Rangeedges));
    else
        histmaton=(histcounts2(mouse_spd(dark&state&laseron), range(dark&state&laseron), Speededges, Rangeedges));
        histmatoff=(histcounts2(mouse_spd(dark&state&~laseron), range(dark&state&~laseron), Speededges, Rangeedges));
    end
    pcolor(histx(Rangeedges),Speededges,(histmaton-histmatoff));
    shading interp
    xlabel('Range, cm')
    ylabel('Speed, cm/s')
    title(sprintf('%s laser on-off (dark)', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'range');
    y=getfield(thresholds, statename{:}, 'mouse_spd');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');
    cl=clim;
    clim([-1 1]*max(abs(cl)))

    nexttile
    if dosqrt
        histmaton=sqrt(histcounts2(mouse_spd(dark&state&laseron), az(dark&state&laseron), Speededges, Azimuthedges));
        histmatoff=sqrt(histcounts2(mouse_spd(dark&state&~laseron), az(dark&state&~laseron), Speededges, Azimuthedges));
    else
        histmaton=(histcounts2(mouse_spd(dark&state&laseron), az(dark&state&laseron), Speededges, Azimuthedges));
        histmatoff=(histcounts2(mouse_spd(dark&state&~laseron), az(dark&state&~laseron), Speededges, Azimuthedges));
    end
    pcolor(histx(Azimuthedges),Speededges,(histmaton-histmatoff));
     colormap blue_to_red
    colormap redbluecmap
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Speed, cm/s')
    title(sprintf('%s laser on-off (dark)', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'az');
    y=getfield(thresholds, statename{:}, 'mouse_spd');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');
    h=colorbar;
    ylabel(h, sprintf('sqrt=%d', dosqrt));
    cl=clim;
    clim([-1 1]*max(abs(cl)))

end

%%%%%%%%%%%%
%%   dark-light
%
% loop through all states

statenames={'hotpursuit', 'chase', 'follow', 'stalk', 'wander', 'pause'};
for statename=statenames

    state=eval(statename{:});
    dosqrt=1;

    figure
    tiledlayout(1,3, "TileSpacing","compact")
    set(gcf, "Position", [440 880 1300 420])
    nexttile
    if dosqrt
        histmaton=sqrt(histcounts2(az(dark&state), range(dark&state), Azimuthedges, Rangeedges));
        histmatoff=sqrt(histcounts2(az(light&state), range(light&state), Azimuthedges, Rangeedges));
    else
        histmaton=(histcounts2(az(dark&state), range(dark&state), Azimuthedges, Rangeedges));
        histmatoff=(histcounts2(az(light&state), range(light&state), Azimuthedges, Rangeedges));
    end
    pcolor(histx(Azimuthedges),histx(Rangeedges),(histmaton'-histmatoff'));
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Range, cm')
    title(sprintf('%s dark-light', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'az');
    y=getfield(thresholds, statename{:}, 'range');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');
    cl=clim;
    clim([-1 1]*max(abs(cl)))

    nexttile
    if dosqrt
        histmaton=sqrt(histcounts2(mouse_spd(dark&state), range(dark&state), Speededges, Rangeedges));
        histmatoff=sqrt(histcounts2(mouse_spd(light&state), range(light&state), Speededges, Rangeedges));
    else
        histmaton=sqrt(histcounts2(mouse_spd(dark&state), range(dark&state), Speededges, Rangeedges));
        histmatoff=sqrt(histcounts2(mouse_spd(light&state), range(light&state), Speededges, Rangeedges));
    end
    pcolor(histx(Rangeedges),histx(Speededges),(histmaton-histmatoff));
    shading interp
    xlabel('Range, cm')
    ylabel('Speed, cm/s')
    title(sprintf('%s dark-light', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'range');
    y=getfield(thresholds, statename{:}, 'mouse_spd');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');
    cl=clim;
    clim([-1 1]*max(abs(cl)))

    nexttile
    if dosqrt
        histmaton=(histcounts2(mouse_spd(dark&state), az(dark&state), Speededges, Azimuthedges));
        histmatoff=(histcounts2(mouse_spd(light&state), az(light&state), Speededges, Azimuthedges));
    else
        histmaton=sqrt(histcounts2(mouse_spd(dark&state), az(dark&state), Speededges, Azimuthedges));
        histmatoff=sqrt(histcounts2(mouse_spd(light&state), az(light&state), Speededges, Azimuthedges));
    end
    pcolor(histx(Azimuthedges),histx(Speededges),(histmaton-histmatoff));
    colormap blue_to_red
    shading interp
    xlabel('Azimuth, degrees')
    ylabel('Speed, cm/s')
    title(sprintf('%s dark-light', statename{:}))
    hold on;
    x=getfield(thresholds, statename{:}, 'az');
    y=getfield(thresholds, statename{:}, 'mouse_spd');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', 'r', 'LineWidth', 1, 'FaceColor', 'none');
    h=colorbar;
    ylabel(h, sprintf('sqrt=%d', dosqrt));
    cl=clim;
    clim([-1 1]*max(abs(cl)))
end


%% test for significant effects of light vs dark on overall geometries
[p,h,stats] =ranksum(  mouse_spd(dark), mouse_spd(light));
fprintf('\nmedian mouse_spd(dark): %.1f, mouse_spd(light): %.1f cms/s', nanmedian(mouse_spd(dark)), nanmedian(mouse_spd(light)))
fprintf('\tranksum p=%.4f z=%.0f', p, stats.zval)

[p,h,stats] =ranksum(  az(dark), az(light));
fprintf('\nmedian az(dark): %.1f, az(light): %.1f cms/s', nanmedian(az(dark)), nanmedian(az(light)))
fprintf('\tranksum p=%.4f z=%.0f', p, stats.zval)

[p,h,stats] =ranksum(  range(dark), range(light));
fprintf('\nmedian range(dark): %.1f, range(light): %.1f cms/s', nanmedian(range(dark)), nanmedian(range(light)))
fprintf('\tranksum p=%.4f z=%.0f', p, stats.zval)

%Molly points out that using frames as the n here is totally invalid.
%if you use mouse as the n, the signranks are not significant.
%if you use session as the n, the p-values are tiny, but is that a valid n?

figure
    tiledlayout(1,3, "TileSpacing","compact")
    set(gcf, "Position", [440 880 1300 420])
nexttile
gold=[    0.9290    0.6940    0.1250];
ndark=histcounts(mouse_spd(dark), Speededges);
nlight=histcounts(mouse_spd(light), Speededges);
p=plot(histx(Speededges), ndark, 'k', histx(Speededges), nlight, 'y');
p(2).Color=gold;
xlabel('speed')

nexttile
ndark=histcounts(az(dark), Azimuthedges, 'norm', 'pdf');
nlight=histcounts(az(light), Azimuthedges, 'norm', 'pdf');
p=plot(histx(Azimuthedges), ndark, 'k', histx(Azimuthedges), nlight, 'y');
p(2).Color=gold;
xlabel('az')

nexttile
ndark=histcounts(range(dark), Rangeedges);
nlight=histcounts(range(light), Rangeedges);
p=plot(histx(Rangeedges), ndark, 'k', histx(Rangeedges), nlight, 'y');
p(2).Color=gold;
xlabel('range')

%%%%%%%%


fprintf('\n there were %d pauses in the dark, for a total pause duration of %d frames', sum(dark(pause_start_frames)), sum(dark&pause))
fprintf('\n there were %d pauses in the light, for a total pause duration of %d frames', sum(light(pause_start_frames)), sum(light&pause))
fprintf('\nthat''s %.1fx more pauses in the dark, lasting %.1fx as long', sum(dark(pause_start_frames))/sum(light(pause_start_frames)), sum(dark&pause)/sum(light&pause))
fprintf('\n(%.1fx more pauses in the light, lasting %.1fx as long)', sum(light(pause_start_frames))/sum(dark(pause_start_frames)), sum(light&pause)/sum(dark&pause))
fprintf('\n')

fprintf('\n there were %d stalks in the dark, for a total stalk duration of %d frames', sum(dark(stalk_start_frames)), sum(dark&stalk))
fprintf('\n there were %d stalks in the light, for a total stalk duration of %d frames', sum(light(stalk_start_frames)), sum(light&stalk))
fprintf('\nthat''s %.1fx more stalks in the dark, lasting %.1fx as long', sum(dark(stalk_start_frames))/sum(light(stalk_start_frames)), sum(dark&stalk)/sum(light&stalk))
fprintf('\n(%.1fx more stalks in the light, lasting %.1fx as long)', sum(light(stalk_start_frames))/sum(dark(stalk_start_frames)), sum(light&stalk)/sum(dark&stalk))

%% plot a single trial
close all
figure
set(gcf, 'pos', [613         818        1796         420])
for i=212 % 144 between 1 and height(metadata)

    fname=metadata{i, 'filename'};
    trialframes=find(contains(filename, fname));
    cricketdrop = metadata{i, 'cricketdrop'};
    captureframe = metadata{i, 'captureframe'};
    if  cricketdrop & captureframe
        trialframes=trialframes(cricketdrop:captureframe);
        %frames that are after cricket drop and before captureframe and are on this trial
        if length(trialframes)>30*200

            aztrial=az(trialframes);
            mouse_spdtrial=mouse_spd(trialframes);
            rangetrial=range(trialframes);
            cricket_spdtrial=cricket_spd(trialframes);
            t=1:length(trialframes);
            t=t/200;

            % aztrial=aztrial/max(aztrial);
            % mouse_spdtrial=mouse_spdtrial/max(mouse_spdtrial);
            % rangetrial=3*rangetrial/max(rangetrial);
            % cricket_spdtrial=cricket_spdtrial/max(cricket_spdtrial);

            aztrial=medfilt1(aztrial, 10);
            mouse_spdtrial=medfilt1(mouse_spdtrial, 10);
            rangetrial=medfilt1(rangetrial, 10);
            cricket_spdtrial=medfilt1(cricket_spdtrial, 10);


            %p=plot(t, mouse_spdtrial, t, cricket_spdtrial, t, rangetrial, t, aztrial);
            %legend('mouse speed', 'cricket speed', 'range', 'az')
            gold=[    0.9290    0.6940    0.1250];
            purple= [0.4940 0.1840 0.5560];

            yyaxis left
            set(gca, 'YColor', gold)
            p1=plot(  t,  rangetrial);
            set(p1, 'linew', 2, 'linestyle', '-')
            set(p1, 'Color', gold)
            ylabel('Range, cm')
            yyaxis right
            p2=plot( t, cricket_spdtrial,t,  aztrial);
            set(p2, 'linew', 2)
            set(p2, 'linew', 2, 'linestyle', '-')
            set(p2(2), 'linew', .5, 'linestyle', '-', 'Color', purple)

            legend('range', 'cricket speed', 'az')
            ylim([0 300])
            yticks([0 90 180])
            ylabel('Azimuth/Cricket Speed')
            title(i)
            set(gca, 'fontsize', 18)
            builtin('pause', 2)
        end
    end
end