function [out,trace]=v3Replay(s,eph,c,independent,scope)
%V3REPLAY Fixed impulses, no shooting-node resets, optional independent targets.
if nargin<4, independent=false; end
if nargin<5, scope='full'; end
assert(ismember(scope,{'full','prefix_witnesses'}),'ctocscreen:v3:replayScope','Unknown replay scope.');
prefix=strcmp(scope,'prefix_witnesses');
out=struct('passed',false,'status','invalid','total_dv_km_s',Inf, ...
 'visit_count',0,'distance_km',inf(35,1),'witness_times_s',nan(35,1), ...
 'height_passed',false,'min_altitude_lower_km',Inf,'sampled_min_altitude_km',Inf, ...
 'validation_level','proposal','failure_reason','');
trace=struct('arcs',{{}},'times',[],'states',[],'heights',{{}},'encounters',{{}});
out.scope=scope; out.targets_independently_propagated=independent&&~prefix;
out.dataset_kind='competition_targets_nominal_model';
if isfield(eph,'synthetic')&&eph.synthetic, out.dataset_kind='synthetic_test_only'; end
try
 m=eph.model; s=ctocscreen.v3Normalize(s,m); x=ctocscreen.initialState(s.initial_q,m.mu,m.re).';
 tq=@(ids,t,mode)ctocscreen.v3QueryTargets(eph,ids,t,mode);
 if independent&&~prefix
  opt=odeset('RelTol',c.verify_reltol,'AbsTol',repmat([1e-12;1e-12;1e-12;1e-15;1e-15;1e-15],35,1), ...
   'MaxStep',min(60,c.max_step_s));
  targetSol=ode89(@targetRhs,[0 s.duration_s],reshape(eph.states0.',[],1),opt);
  assert(targetSol.x(end)>=s.duration_s,'ctocscreen:v3:targetIntegration','Target integration stopped.');
  tq=@independentQuery;
 end
 knots=unique([0;s.maneuver_times_s;s.duration_s]); height=true;
 trace.times=knots; trace.states=zeros(numel(knots),6);
 for k=1:numel(knots)
  b=find(s.maneuver_times_s==knots(k));
  if ~isempty(b), x(4:6)=x(4:6)+s.delta_v_km_s(b,:).'; end
  trace.states(k,:)=x.';
  if k==numel(knots), break; end
  [x,~,sol]=ctocscreen.v3Arc(x,knots(k),knots(k+1),m,c,false,independent);
  trace.arcs{end+1}=sol;
  h=ctocscreen.v3Height(sol,m,c); height=height&&h.passed;
  trace.heights{end+1}=h;
  out.min_altitude_lower_km=min(out.min_altitude_lower_km,h.lower_bound_altitude_km);
  out.sampled_min_altitude_km=min(out.sampled_min_altitude_km,h.sampled_min_altitude_km);
  if prefix
   % Construction already found these witnesses. Recheck them without another
   % all-target closest-approach search; this mode cannot certify the mission.
   z=struct('distance_km',inf(35,1),'time_s',nan(35,1));
   ids=find(isfinite(s.witness_times_s)&s.witness_times_s>=sol.x(1)&s.witness_times_s<=sol.x(end));
   if ~isempty(ids)
    times=s.witness_times_s(ids); y=deval(sol,times); r=tq(ids,times,'pairs');
    z.distance_km(ids)=vecnorm(y(1:3,:).'-r,2,2); z.time_s(ids)=times;
   end
  else
   z=ctocscreen.v3ScanArc(sol,tq,c,s.witness_times_s);
  end
  trace.encounters{end+1}=z;
  better=z.distance_km<out.distance_km;
  out.distance_km(better)=z.distance_km(better); out.witness_times_s(better)=z.time_s(better);
 end
 out.final_state=x.'; out.height_passed=height;
 out.total_dv_km_s=sum(vecnorm(s.delta_v_km_s,2,2)); out.visit_count=sum(out.distance_km<=1);
 out.inclination_penalty=ctocscreen.v3TracePlanePenalty(s,trace,c);
 out.passed=height&&out.visit_count==35; out.status='constraints_not_met';
 out.validation_level='nominal_j2_screened';
 if out.passed, out.status='complete'; end
 if independent, out.validation_level='nominal_j2_independent_checked'; end
 if independent&&out.passed, out.validation_level='nominal_j2_independent_verified'; end
 if prefix, out.validation_level='nominal_j2_prefix_checked'; end
 out.official_alignment_verified=m.official_alignment_verified;
 out.event_scope='positive witnesses for all targets; no claim of exhaustive absence';
 if prefix, out.event_scope='supplied witnesses only; target ephemeris reused; not full independent verification'; end
catch err
 out.status='propagation_failure'; out.failure_reason=err.message; out.failure_id=err.identifier;
 out.passed=false; out.total_dv_km_s=Inf;
end
 function dx=targetRhs(t,y)
  z=reshape(y,6,35).'; a=ctocscreen.v3ReferenceForce(t,z(:,1:3),m);
  dx=reshape([z(:,4:6),a].',[],1);
 end
 function [r,v]=independentQuery(ids,t,mode)
  if isempty(ids), ids=1:35; end
  ids=ids(:); t=t(:); ni=numel(ids); nt=numel(t);
  if strcmp(mode,'grid')
   ids=repmat(ids,nt,1); t=repelem(t,ni);
  elseif isscalar(ids)
   ids=repmat(ids,nt,1);
  elseif isscalar(t)
   t=repmat(t,ni,1);
  end
  [ut,~,it]=unique(t); y=deval(targetSol,ut); r=zeros(numel(ids),3); v=r;
  for ii=1:numel(ids)
   row=6*(ids(ii)-1); r(ii,:)=y(row+(1:3),it(ii)).'; v(ii,:)=y(row+(4:6),it(ii)).';
  end
  if strcmp(mode,'grid')
   r=permute(reshape(r,ni,nt,3),[1 3 2]); v=permute(reshape(v,ni,nt,3),[1 3 2]);
  end
 end
end
