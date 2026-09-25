function cache = v3LoadTargetEphemeris(filename,requireOfficialAlignment)
%V3LOADTARGETEPHEMERIS Load once per worker and reject stale input/code caches.
% By default allows explicitly labelled nominal orientation for local work.
if nargin<2, requireOfficialAlignment=false; end
filename=char(filename);
assert(isfile([filename '.sha256']),'ctocscreen:v3:checksum','Missing cache checksum.');
assert(strcmp(strtrim(fileread([filename '.sha256'])),ctocscreen.v3FileHash(filename)), ...
 'ctocscreen:v3:checksum','Cache bytes changed or incomplete.');
s=load(filename,'cache'); cache=s.cache;
assert(strcmp(cache.schema_version,'target_ephemeris_v3_1'),'ctocscreen:v3:schema','Wrong schema.');
assert(isequaln(cache.signature,ctocscreen.v3PreprocessSignature()), ...
 'ctocscreen:v3:staleCache','Source/data/version changed. Rebuild preprocessing cache.');
assert(cache.validation.passed,'ctocscreen:v3:unvalidated','Numerical validation failed.');
assert(~requireOfficialAlignment || cache.model.official_alignment_verified, ...
 'ctocscreen:v3:alignmentPending','ATK native frame/EOP alignment is not certified.');
end
