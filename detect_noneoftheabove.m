function [noneoftheabove, noneoftheabove_start_frames, noneoftheabove_end_frames]=detect_noneoftheabove(varargin) 

% usage: [noneoftheabove, noneoftheabove_start_frames, noneoftheabove_end_frames]=detect_noneoftheabove(state1, state2, ...) 

% detects noneoftheabove state,  a continuous binary state variable,
% defined by the absence of any of the states passed as input
% input states should be logicals of length 1:numframes



fprintf('\n\ndetecting "none of the above" state... ')
if nargin==0 return, end

states=zeros(nargin, length(varargin{1}));
for n=1:nargin
    states(n,:)=varargin{n};
end

noneoftheabove=~sum(states);
noneoftheabove=noneoftheabove(:);

noneoftheabove_start_frames=find(diff(noneoftheabove)==1);
noneoftheabove_end_frames=find(diff(noneoftheabove)==-1);
if noneoftheabove(end)==1
    noneoftheabove_end_frames=[noneoftheabove_end_frames; length(noneoftheabove)];
end
if noneoftheabove(1)==1
    noneoftheabove_start_frames=[1; noneoftheabove_start_frames];
end

fprintf('\nfound %d noneoftheabove epochs, total of %d frames (%.0f%%)', length(noneoftheabove_start_frames), sum(noneoftheabove), 100*sum(noneoftheabove)/length(noneoftheabove))

