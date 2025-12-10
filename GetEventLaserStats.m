function [event_counts_ON, event_counts_OFF, on_pvalue, RateRatio_ON_vs_OFF] = GetEventLaserStats(event_frames, event_type, metadata, filename, verbose)
% compute laser on/off event rate statistics using a Poisson Regression Model
% usage: [event_counts_ON, event_counts_OFF, on_pvalue, RateRatio_ON_vs_OFF] = GetEventLaserStats(event_frames, event_type, metadata, filename, verbose);
% set 'verbose' to 1 to print all the statistics output, or 0 to run
% silently and just return sign and p-value

%to do: write a version of this that separates light/dark trials

%use metadata to compile event and trial counts with/without laser
event_counts_OFF = 0;
event_counts_ON = 0;
for eventframe=event_frames(:)'; %at 200 fps
    laserval = metadata{contains(metadata.filename, filename{eventframe}), 'laser_value'};
    if ~isnan(laserval) %no idea why there are a couple nan laser values
        if laserval
            event_counts_ON=event_counts_ON+1;
        else
            event_counts_OFF=event_counts_OFF+1;
        end
    end
end

% note that it would be wrong to count number of ON or OFF trials or frames in the above loop,
% because the loop skips the trials on which there aren't any events
% instead we just count all laser on and off trials from the metadata, below (it's fast)
laservals=metadata.laser_value;
laservals=laservals(~isnan(laservals)); %remove nans
trialnums_ON=find(laservals);
trialnums_OFF=find(~laservals);
numtrials_ON=sum(laservals);
numtrials_OFF=sum(~laservals);
framecounts_ON=[];
framecounts_OFF=[];
for trial=trialnums_ON
    numframes_in_trial=metadata{trial, 'captureframe'} - metadata{trial, 'cricketdrop'};
    framecounts_ON = [framecounts_ON numframes_in_trial];
end
for trial=trialnums_OFF
    numframes_in_trial=metadata{trial, 'captureframe'} - metadata{trial, 'cricketdrop'};
    framecounts_OFF = [framecounts_OFF numframes_in_trial];
end
totalnumframesON=sum(framecounts_ON);
totalnumframesOFF=sum(framecounts_OFF);


if verbose
    fprintf('\n%d ON trials, for a total of %d ON frames', numtrials_ON, totalnumframesON);
    fprintf('\n%d OFF trials, for a total of %d OFF frames', numtrials_OFF, totalnumframesOFF);

    %ways we could normalize: by numtrials, by numframes, or by events/minute

    fprintf('\n%s, laser ON: %d events, laser OFF: %d events ', event_type, event_counts_ON, event_counts_OFF);
    fprintf('\n%s, laser ON: %.2f events/trial, laser OFF: %.2f events/trial ', event_type, event_counts_ON/numtrials_ON, event_counts_OFF/numtrials_OFF);
    fprintf('\n%s, laser ON: %.2e events/frame, laser OFF: %.2e events/frame ', event_type, event_counts_ON/totalnumframesON, event_counts_OFF/totalnumframesOFF);
    fprintf('\n%s, laser ON: %.2f events/minute, laser OFF: %.2f events/minute ', event_type, event_counts_ON*200*60/totalnumframesON, event_counts_OFF*200*60/totalnumframesOFF);
end


%% Poisson Regression Model

% 1. Define the number of trials and events for each condition.
Trials = [numtrials_ON; numtrials_OFF];
Events = [event_counts_ON; event_counts_OFF];

% 2. Define the condition labels (categorical is best for GLM)
Condition = categorical({'ON'; 'OFF'});

% 3. Create a table for cleaner data management
T = table(Condition, Trials, Events);

% Display the simple data structure
if verbose
    fprintf('\nPoisson Regression test for effect of laser on rate of %s events:\n', event_type);
    disp('--- Input Data ---');
    disp(T);
end

% Poisson Regression Model

% The Generalized Linear Model (GLM) formula is:
% log(Expected_Events) = beta0 + beta1 * Condition(ON vs OFF) + log(Trials)
% The log(Trials) term is the OFFSET, which converts the model from
% predicting the COUNT to predicting the RATE (Count/Trial).

% Fit the Poisson GLM model
model = fitglm(T, 'Events ~ Condition', ...
    'Distribution', 'poisson', ...
    'Offset', log(T.Trials));

% Interpretation

% The coefficient for Condition_OFF represents the log-rate difference
% between the OFF condition and the reference ON condition.

if verbose
    disp('--- GLM Results: Comparing Log-Rates ---');
    disp(model);
end

% RateRatio_ON_vs_OFF (e^beta) tells you how much more likely an event is
% in the ON condition relative to the OFF condition.

% Find the coefficient for the Condition (e.g., the ON vs OFF effect)
coeff_table = model.Coefficients;
on_coeff = coeff_table{'Condition_ON', 'Estimate'};
on_pvalue = coeff_table{'Condition_ON', 'pValue'};

RateRatio_ON_vs_OFF = exp(on_coeff);

if verbose
    fprintf('The rate ratio (ON vs OFF) is: %.3f, ', RateRatio_ON_vs_OFF);

    if on_pvalue < 0.05
        fprintf('which is statistically significant (p = %.4f).\n', on_pvalue);
    else
        fprintf('which is not statistically significant (p = %.4f).\n', on_pvalue);
    end
end
