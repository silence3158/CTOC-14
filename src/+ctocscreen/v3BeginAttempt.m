function [state,source]=v3BeginAttempt(state,hasWork,hasRecovery,hasInitial)
%V3BEGINATTEMPT Reserve a source slot before any operation can fail or skip.
source=ctocscreen.v3ChooseSource(numel(state.source_history),hasWork,hasRecovery);
if hasInitial, source='initial'; end
ordinal=1+sum(cellfun(@(h)strcmp(h.source,source),state.source_history));
state.source_history{end+1}=struct('source',source,'outcome','started','source_attempt',ordinal);
end
