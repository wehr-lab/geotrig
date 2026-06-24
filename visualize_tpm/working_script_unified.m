%  working_script_unified

% Usage:
%
%
% plot_unified_interaction(TPM, statenames, hawkesParams, eventnames, ...
%     stateRates, glmResults, eventStateJS, ...
%     'TPM_thresh', 0.08, ...
%     'Hawkes_thresh', 0.08, ...
%     'StateEvent_thresh', 0.8, ...
%     'EventState_thresh', 0.08);

% hawkes_params are loaded from geotrig/datacache
  load /Users/wehr/Documents/Analysis/geotrig/datacache 

TPM=TPMdark;
statenames = {'hot pursuit','chase','following','stalk','wander','pause'};
eventnames = {'failed_approach','contact_loss', ...
    'contact_gain','intercept','cricket_jump','rangemin'};

nStates = numel(statenames);
nEvents = numel(eventnames);

eventFrames={failed_approach_event_frames,contact_loss_event_frames, ...
    contact_gain_event_frames,intercept_event_frames,cricket_jump_event_frames,rangemin_event_frames};
stateMask = [    hotpursuit(:),chase(:),follow(:),...
    stalk(:), wander(:), pause(:)];
fps=200;

if 0
    %% generate the stuff

    [stateRates, rateMat] = pp_state_rates(eventFrames, eventnames, stateMask, statenames, fps);

    glmResults = pp_glm(eventFrames, eventnames, stateMask, statenames, num_frames, fps);

    Event_conditioned_TPMs;% run to return results struct
    %Assemble eventStateJS from your existing results struct:
    eventStateJS = zeros(nEvents, nStates);
    for e = 1:nEvents
        eventStateJS(e,:) = results(e).jsdiv;
    end
end

%% 
plot_unified_interaction(TPM, statenames, hawkes_params, eventnames, ...
    stateRates, glmResults, eventStateJS, ...
    'TPM_thresh', 0.08, ... %black
    'Hawkes_thresh', 0.08, ... %orange
    'StateEvent_thresh', 1, ... %blue
    'EventState_thresh', 0.02); %red




