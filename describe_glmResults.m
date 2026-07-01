% convert GLM results into plain english


if ~exist('glmResults')
    [glmResults] = pp_glm(eventFrames, eventnames, stateMask, statenames, ...
    num_frames, fps);
end

fprintf('\nsignificant coeeficients between event and event/lag/state, ranked by p-value (<.00001)')
for e=1:6

 fprintf('\n___________________________')
 fprintf('\n%s:',    glmResults(e).event   )
 fprintf('\n___________________________')
    [pvals, order]=sort(glmResults(e).pvals);

    coefs=glmResults(e).coefs(order);
    sigidx=find(pvals<.00001);
    featNames=glmResults(e).featNames(order);

    for c=sigidx'
        fprintf('\n%s:\n\tp=%.4f, \tcoef=%.2f', featNames{c}, pvals(c), coefs(c))
    end
end