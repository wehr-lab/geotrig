% draw_thresholds
% draw boxes to show state thresholds
% and think about how to resolve potential overlaps

statenames={'chase1','chase2','chase3','chase4', 'wander', 'pause', 'stalk'};

thresholds.chase.az=[0 30];
thresholds.chase.mouse_spd=[10 60];
thresholds.chase.range=[5 10];

thresholds.chase1.mouse_spd=[5 10];
thresholds.chase2.mouse_spd=[10 20];
thresholds.chase3.mouse_spd=[20 30];
thresholds.chase4.mouse_spd=[30 50];
%speedvalues = [ 5 10 20 30 50];

%copy az and range from chase
thresholds.chase1.az=thresholds.chase.az;
thresholds.chase2.az=thresholds.chase.az;
thresholds.chase3.az=thresholds.chase.az;
thresholds.chase4.az=thresholds.chase.az;
thresholds.chase1.range=thresholds.chase.range;
thresholds.chase2.range=thresholds.chase.range;
thresholds.chase3.range=thresholds.chase.range;
thresholds.chase4.range=thresholds.chase.range;


thresholds.wander.az=[30 180];
thresholds.wander.mouse_spd=[5 15];
thresholds.wander.range=[15 60];

thresholds.pause.az=[0 180]; %wide open
thresholds.pause.mouse_spd=[0 5];
thresholds.pause.range=[0 60]; %wide open

thresholds.stalk.az=[30 180]; %wide open
thresholds.stalk.mouse_spd=[3 15];
thresholds.stalk.range=[5 40]; %wide open


Azimuthedges=linspace(0, 180, 50);
Rangeedges=linspace(0, 60, 50);
Speededges=linspace(0, 80, 50);


figure
tiledlayout(1,3, "TileSpacing","compact")
set(gcf, "Position", [440 880 1300 1*420])

cmap=prism(length(statenames));

nexttile
xlabel('Azimuth, degrees')
ylabel('Range, cm/s')
hold on;
i=0;
for statename=statenames
i=i+1;
    %state=eval(statename{:});

    x=getfield(thresholds, statename{:}, 'az');
    y=getfield(thresholds, statename{:}, 'range');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', cmap(i,:), 'LineWidth', 1, 'FaceColor', 'none');
    t=text(mean(x), mean(y), statename, 'color', cmap(i,:));
end

nexttile
xlabel('Range, degrees')
ylabel('Speed, cm/s')
hold on;
i=0;
for statename=statenames
i=i+1;
    %state=eval(statename{:});

    x=getfield(thresholds, statename{:}, 'range');
    y=getfield(thresholds, statename{:}, 'mouse_spd');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', cmap(i,:), 'LineWidth', 1, 'FaceColor', 'none');
    t=text(mean(x), mean(y), statename, 'color', cmap(i,:));
end

nexttile
xlabel('Azimuth, degrees')
ylabel('Speed, cm/s')
hold on;
i=0;
for statename=statenames
i=i+1;
    %state=eval(statename{:});

    x=getfield(thresholds, statename{:}, 'az');
    y=getfield(thresholds, statename{:}, 'mouse_spd');
    pos = [x(1), y(1), x(2)-x(1), y(2)-y(1)];
    rectangle('Position', pos, 'EdgeColor', cmap(i,:), 'LineWidth', 1, 'FaceColor', 'none');
    t=text(mean(x), mean(y), statename, 'color', cmap(i,:));
end



