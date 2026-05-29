function [event_counts_1, event_counts_2, on_pvalue, RateRatio_1_vs_2, totalnumframes1, totalnumframes2] = GetEventLaserStats_kip02(event_frames, event_type, metadata, filename, condition1, condition2, verbose)
% compute laser on/off event rate statistics using a Poisson Regression Model
% usage: [event_counts_ON, event_counts_OFF, on_pvalue, RateRatio_ON_vs_OFF, totalnumframesON, totalnumframesOFF] = GetEventLaserStats(event_frames, event_type, metadata, filename, condition1, condition2, verbose);
% condition1, condition2 are the 5-digit conditions to be compared
% set 'verbose' to 1 to print all the statistics output, or 0 to run
% silently and just return sign and p-value

%to do: write a version of this that separates light/dark trials

%use metadata to compile event and trial counts with condition1 or condition2
event_counts_1 = 0;
event_counts_2 = 0;
for eventframe=event_frames(:)'     % at 200 fps
    condition = metadata{contains(metadata.filename, filename{eventframe}), 'condition'};
    if ~isnan(condition)
        if condition==condition1
            event_counts_1 =event_counts_1 +1;
        elseif condition==condition2
            event_counts_2=event_counts_2+1;
        end
    end
end

% note that it would be wrong to count number of ON or OFF trials or frames in the above loop,
% because the loop skips the trials on which there aren't any events
% instead we just count all laser on and off trials from the metadata, below (it's fast)
conditions=metadata.condition;
conditions=conditions(~isnan(conditions)); %remove nans
trialnums_1=find(conditions == condition1);
trialnums_2=find(conditions == condition2);
numtrials_1=sum(trialnums_1);
numtrials_2=sum(trialnums_2);
framecounts_1=[];
framecounts_2=[];
for trial=trialnums_1
    numframes_in_trial=metadata{trial, 'captureframe'} - metadata{trial, 'cricketdrop'};
    framecounts_1 = [framecounts_1 numframes_in_trial];
end
for trial=trialnums_2
    numframes_in_trial=metadata{trial, 'captureframe'} - metadata{trial, 'cricketdrop'};
    framecounts_2 = [framecounts_2 numframes_in_trial];
end
totalnumframes1=sum(framecounts_1);
totalnumframes2=sum(framecounts_2);


if verbose
    fprintf('\n%d condition1 trials, for a total of %d condition1 frames', numtrials_1, totalnumframes1);
    fprintf('\n%d condition2 trials, for a total of %d condition2 frames', numtrials_2, totalnumframes2);

    %ways we could normalize: by numtrials, by numframes, or by events/minute

    fprintf('\n%s, condition1: %d events, condition2: %d events ', event_type, event_counts_1, event_counts_2);
    fprintf('\n%s, condition1: %.2f events/trial, condition2: %.2f events/trial ', event_type, event_counts_1/numtrials_1, event_counts_2/numtrials_2);
    fprintf('\n%s, condition1: %.2e events/frame, condition2: %.2e events/frame ', event_type, event_counts_1/totalnumframes1, event_counts_2/totalnumframes2);
    fprintf('\n%s, condition1: %.2f events/minute, condition2: %.2f events/minute ', event_type, event_counts_1*200*60/totalnumframes1, event_counts_2*200*60/totalnumframes2);
end


%% Poisson Regression Model

% 1. Define the number of trials and events for each condition.
% altered to use the #frames for each condition

fprintf('modified here to use #frames/condition in the model\n')

Frames = [totalnumframes1; totalnumframes2];
%Trials = [numtrials_ON; numtrials_OFF];
Events = [event_counts_1; event_counts_2];

% 2. Define the condition labels (categorical is best for GLM)
Condition = categorical({'condition1'; 'condition2'});

% 3. Create a table for cleaner data management
%T = table(Condition, Trials, Events);
T = table(Condition, Frames, Events);

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
% model = fitglm(T, 'Events ~ Condition', ...
%     'Distribution', 'poisson', ...
%     'Offset', log(T.Trials));

model = fitglm(T, 'Events ~ Condition', ...
    'Distribution', 'poisson', ...
    'Offset', log(T.Frames));

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
on_coeff = coeff_table{'Condition_condition2', 'Estimate'};
on_pvalue = coeff_table{'Condition_condition2', 'pValue'};

RateRatio_1_vs_2 = exp(on_coeff);

if verbose
    fprintf('The rate ratio (1 vs 2) is: %.3f, ', RateRatio_1_vs_2);

    if on_pvalue < 0.05
        fprintf('which is statistically significant (p = %.4f).\n', on_pvalue);
    else
        fprintf('which is not statistically significant (p = %.4f).\n', on_pvalue);
    end
end
