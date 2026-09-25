function s=v3Normalize(s,m)
%V3NORMALIZE Validate free schedule; exact simultaneous impulses are summed.
assert(strcmp(s.schema_version,'free_maneuver_v3') && strcmp(s.dynamics_id,'central_j2'), ...
 'ctocscreen:v3:schedule','Wrong schedule interface.');
q=s.initial_q(:).'; validateattributes(q,{'double'},{'numel',6,'finite','real'});
assert(q(1)>=m.re+590 && q(1)<=m.re+610 && hypot(q(2),q(3))<.001 && q(4)>=0 && q(4)<=pi, ...
 'ctocscreen:v3:initialOrbit','Initial orbit outside formal bounds.');
s.initial_q=q; t=s.maneuver_times_s(:); d=s.delta_v_km_s;
if ~isfield(s,'root_key'), s.root_key=sprintf('%.17g/',q); end
assert(isscalar(s.duration_s)&&isfinite(s.duration_s)&&s.duration_s>0&&s.duration_s<=m.horizon_s);
assert(size(d,1)==numel(t)&&size(d,2)==3&&all(isfinite(d(:)))&&all(isfinite(t)));
assert(all(t>=0 & t<=s.duration_s),'ctocscreen:v3:schedule','Burn outside task.');
[t,~,group]=unique(t); merged=zeros(numel(t),3);
for k=1:numel(group), merged(group(k),:)=merged(group(k),:)+d(k,:); end
keep=any(merged~=0,2); s.maneuver_times_s=t(keep); s.delta_v_km_s=merged(keep,:);
assert(numel(s.witness_times_s)==35,'ctocscreen:v3:schedule','Need 35 witness slots (NaN for unknown).');
s.witness_times_s=s.witness_times_s(:);
w=s.witness_times_s;
assert(all(isnan(w)|(isfinite(w)&w>=0&w<=s.duration_s)),'ctocscreen:v3:schedule','Invalid witness times.');
if isfield(s,'visit_plan_times_s')
 w=s.visit_plan_times_s(:);
 assert(numel(w)==35&&all(isnan(w)|(isfinite(w)&w>=0&w<=s.duration_s)), ...
  'ctocscreen:v3:schedule','Invalid proposed visit times.');
 s.visit_plan_times_s=w;
end
end
