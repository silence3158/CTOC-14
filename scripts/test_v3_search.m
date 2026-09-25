function results=test_v3_search()
%TEST_V3_SEARCH Synthetic correctness tests; no competition performance claim.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'),fullfile(sim,'tests'));
suite=testsuite({fullfile(sim,'tests','TestV3Search.m'),fullfile(sim,'tests','TestV3RecoveryFlow.m'), ...
 fullfile(sim,'tests','TestV3Diversity.m'),fullfile(sim,'tests','TestV3Construction.m'), ...
 fullfile(sim,'tests','TestV3Replan.m'),fullfile(sim,'tests','TestV3ExecutionGates.m'), ...
 fullfile(sim,'tests','TestV3PrefixConsistency.m'),fullfile(sim,'tests','TestV3ColdStart.m')});
results=run(suite);
folder=fullfile(sim,'runs','v3','development'); if ~isfolder(folder), mkdir(folder); end
save(fullfile(folder,'test_results.mat'),'results');
assert(all([results.Passed]),'ctocscreen:v3:tests','V3 correctness tests failed.');
end
