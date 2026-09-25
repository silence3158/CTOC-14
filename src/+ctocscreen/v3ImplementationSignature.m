function s=v3ImplementationSignature(eph)
%V3IMPLEMENTATIONSIGNATURE Versioned search checkpoint dependencies.
folder=fileparts(mfilename('fullpath')); files=dir(fullfile(folder,'v3*.m'));
extra={'initialState.m','enumerateBranches.m','stumpff.m','propagateTwoBody.m','checkArc.m', ...
 'private/lambert_izzo_gooding.m'};
names=[{files.name},extra]; s.names=names; s.sha256=cell(size(names));
for k=1:numel(names), s.sha256{k}=ctocscreen.v3FileHash(fullfile(folder,names{k})); end
s.ephemeris_signature=eph.signature; s.matlab_version=version;
s.schema='free_maneuver_search_v3_26';
end
