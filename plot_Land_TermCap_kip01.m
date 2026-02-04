% quick script to compare cricket land, termCap measured vs calculated
% should start here: X:\PreyCapture\A1Suppression

% also checks and rewrites Segmentation file if incomplete!!!!

% need to have "metadata"
metadataFN = 'metadata_condition_10000p0';
metadata=readtable(metadataFN,'Delimiter',',','ReadVariableNames',true, 'ReadRowNames',false);

% move to data directory
cd save

F2 = figure; hold on
F3 = figure; hold on
BadCricketDrops = 0;
BadCaptures = 0;
for i = 1:height(metadata)
    fprintf('working on %d of %d\n',i,height(metadata))
    cd(metadata.folder_path{i})
    load Segmentation

    if ~isfield(Seg,'FatalFlag')
        [Seg] = checkSegFieldNames(Seg);
        save('Segmentation.mat','Seg')
        fprintf('      Wrote Segmentation file in %s\n',metadata.folder_path{i})
    end

    if Seg.FatalFlag==0
        if ~isempty(Seg.Land)
            BadCricketDrops = BadCricketDrops+1;
            figure(F2)
            plot(metadata.cricketdrop(i),Seg.Land,'b.','markersize',24)
        else
            fprintf('no Landing in %s\n',metadata.folder_path{i})
        end
        if ~isempty(Seg.TerminalCap)
            BadCaptures= BadCaptures+1;
            figure(F3)
            plot(metadata.captureframe(i),Seg.TerminalCap,'r.','markersize',24)
        else
            fprintf('no TerminalCap in %s\n',metadata.folder_path{i})
        end
    end
    clear Seg
    cd ..
end
figure(F2)
title('Cricket drop')
plot([0 65000],[9 65000],'k','linewidth',2)
xlabel('metadata')
ylabel('measured')
BadCricketDrops

figure(F3)
title('Terminal Capture')
plot([0 90000],[9 90000],'k','linewidth',2)
xlabel('metadata')
ylabel('measured')
BadCaptures