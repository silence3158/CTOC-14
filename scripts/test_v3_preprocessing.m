function results=test_v3_preprocessing(folder)
%TEST_V3_PREPROCESSING Regression checks for an explicitly selected cache.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
old=getenv('CTOC_V3_TEST_CACHE'); cleanup=onCleanup(@()setenv('CTOC_V3_TEST_CACHE',old));
setenv('CTOC_V3_TEST_CACHE',fullfile(folder,'target_ephemeris.mat'));
results=runtests(fullfile(sim,'tests','TestV3Preprocessing.m'));
save(fullfile(folder,'test_results.mat'),'results');
assert(all([results.Passed]),'ctocscreen:v3:tests','Preprocessing tests failed.');
end
