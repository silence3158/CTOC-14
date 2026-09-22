function [pl,core]=structureTemplate(source,first,counts,p,cfg)
%STRUCTURETEMPLATE Replace exactly three consecutive S1 arcs, retain all IDs.
old=ctocscreen.arcPlan(source.schedule,source.evaluation);
assert(all(old.counts(first:first+2)==1)&&sum(counts)==3&&counts(end)>0);
assert(all(counts>=0)&&all(counts==floor(counts)));
if ~isfield(old,'aim_offsets_km'),old.aim_offsets_km=zeros(numel(old.counts),3);end
N=numel(old.counts);ends=cumsum(old.counts);core.ids=old.ids(ends(first):ends(first+2));
core.first_event=ends(first);core.last_event=ends(first+2);core.first_arc=first;
core.last_arc=first+numel(counts)-1;pl=old;
newRef=zeros(numel(counts),3);newAim=newRef;waits=zeros(numel(counts),1);zeroDV=newRef;
visited=0;
for m=1:numel(counts)
 if counts(m)>0
  visited=visited+counts(m);anchor=first+visited-1;
  newRef(m,:)=old.reference_v(anchor,:);newAim(m,:)=old.aim_offsets_km(anchor,:);
  if m==1,waits(m)=old.waits(first);else,waits(m)=old.waits(first+visited-counts(m));end
 else
  assert(m==1,'Current experiment templates place zero-visit burn first.');
  waits(m)=old.waits(first);zeroDV(m,:)=source.schedule.delta_v_km_s(first,:);
 end
end
if counts(1)==0
 previous=0;if first>1,previous=old.times(ends(first-1));end
 available=old.times(ends(first))-previous-old.waits(first);
 waits(2)=min(600,.1*available);
 assert(waits(2)>cfg.minimum_gap_s);
 x=source.evaluation.preburn_states(first,:);x(4:6)=x(4:6)+zeroDV(1,:);
 x=ctocscreen.propagateTwoBody(x,waits(2),p.mu_km3_s2);newRef(2,:)=x(4:6);
end
before=1:first-1;after=first+3:N;
pl.counts=[old.counts(before);counts(:);old.counts(after)];
pl.waits=[old.waits(before);waits;old.waits(after)];
pl.reference_v=[old.reference_v(before,:);newRef;old.reference_v(after,:)];
pl.aim_offsets_km=[old.aim_offsets_km(before,:);newAim;old.aim_offsets_km(after,:)];
pl.zero_delta_v=[zeros(numel(before),3);zeroDV;zeros(numel(after),3)];
% Branch references seed a family; after this first selection, the continuous
% call locks the chosen branch labels. The exact global target order is fixed.
[~,r]=ctocscreen.rebuildArcPlan(pl,p,cfg);
if r.completed,pl.locked_branch_ids=r.branch_ids;end
end
