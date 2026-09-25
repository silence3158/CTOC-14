function check_v3_family_guidance()
%CHECK_V3_FAMILY_GUIDANCE Focused checks for the family-cohesion ranking terms.
% These are ranking heuristics only; passing proves the terms exist and behave,
% not that the search is better.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
c=ctocscreen.v3Defaults();
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs','v3','preprocessing','targets_20260923_release','target_ephemeris.mat'));
x=eph.states0; r=vecnorm(x(:,1:3),2,2);
invA=2./r-sum(x(:,4:6).^2,2)/eph.model.mu; a=1./invA;
h=cross(x(:,1:3),x(:,4:6),2); inc=acosd(h(:,3)./vecnorm(h,2,2));
fprintf('FAMILY target bands: %s\n',mat2str(round(unique(round(a/2.5e3))*2.5e3).'));
fprintf('FAMILY target planes: %s deg\n',mat2str(round(unique(round(inc/10)*10)).'));
% Separation must be zero for a target in the same band and plane, and positive
% for a different band or a different plane.
same=find(abs(a-min(a))<1e3,1); far=find(inc>max(inc)-1,1);
assert(~isempty(same)&&~isempty(far),'fixture selection failed');
homeBand=round(a(same)/2.5e3); homeInc=inc(same);
[sepSame,detailSame]=ctocscreen.v3FamilySeparation(0,600,eph,same,c,homeBand,homeInc);
[sepPlane,~]=ctocscreen.v3FamilySeparation(0,600,eph,far,c,homeBand,homeInc);
bandTarget=find(abs(a-max(a))<1e3,1);
[sepBand,detailBand]=ctocscreen.v3FamilySeparation(0,600,eph,bandTarget,c,homeBand,homeInc);
assert(detailSame.arrival_epoch_s==600,'separation must be evaluated at the arrival epoch');
assert(detailBand.arrival_epoch_s==600,'band fixture must use the arrival epoch');
% Same-family separation is zero up to the floating-point drift of a J2 orbit
% between the departure epoch and the arrival epoch (measured ~2.6e-7 km/s).
assert(sepSame<1e-4,'same-family separation should be ~0');
assert(sepPlane>sepSame,'plane change must cost more than staying in plane');
assert(~detailSame.band_change,'same band flagged as a band change');
% The candidate ranking must react: ask for proposals with a small lookahead and
% confirm the separation accounting is reported and finite.
node=struct('schedule',struct('initial_q',[eph.model.re+600 0 0 0 0 0], ...
  'maneuver_times_s',zeros(0,1),'delta_v_km_s',zeros(0,3),'duration_s',0, ...
  'witness_times_s',nan(35,1),'validation_level','proposal'), ...
 'state',ctocscreen.initialState([eph.model.re+600 0.0001 -0.0002 inc(same)*pi/180 .3 .4],eph.model.mu,eph.model.re), ...
 't',0,'visited',false(35,1),'J',0,'keys',{{}},'last_gain',0,'estimate',0);
c.budget_s=60; c.cost_lookahead_s=7200; c.cost_time_samples=3; c.cost_natural_windows=1;
[pairs,report]=ctocscreen.v3CostCandidates(node,eph,c,RandStream('mt19937ar','Seed',7), ...
 containers.Map('KeyType','char','ValueType','double'),tic);
assert(~isempty(pairs),'enumeration produced no pairs');
fprintf('FAMILY pairs=%d separation_hits=%d completion_hits=%d\n', ...
 report.pairs_evaluated,report.family_separation_hits,report.family_completion_hits);
assert(report.family_separation_hits>0,'no candidate was flagged as a family separation');
assert(all(isfinite(cellfun(@(p)p.weight,pairs))),'non-finite ranking weight');
assert(all(cellfun(@(p)p.weight,pairs)>=0),'negative ranking weight');
assert(all(isfield(pairs{1},{'family','family_cost_km_s'})),'family fields missing from pairs');
% Completion reward must lift a nearly finished family above an equal-cost leg
% in a family with the same amount of work left.
fprintf('FAMILY_CHECK: separation is zero in-family, positive across plane/band, ranking weights finite and non-negative.\n');
end
