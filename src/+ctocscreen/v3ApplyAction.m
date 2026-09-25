function child=v3ApplyAction(node,departure,dv,finish,target,eph,c,arc)
%V3APPLYACTION Propagate waiting and departure from the actual prefix state.
if nargin<8, arc=[]; end
assert(departure>=node.t&&finish>departure&&finish<=eph.model.horizon_s);
s=node.schedule; if isfield(s,'visit_plan_times_s'), s=rmfield(s,'visit_plan_times_s'); end
dv=dv(:);
if any(dv~=0)
 assert(numel(s.maneuver_times_s)<c.max_maneuvers&&all(abs(dv)<=c.pulse_component_bound), ...
  'ctocscreen:v3:actionBounds','Configurable maneuver search bounds exceeded.');
 s.maneuver_times_s(end+1,1)=departure; s.delta_v_km_s(end+1,:)=dv.';
end
s.duration_s=finish; xx=node.state(:); start=node.t; visited=node.visited; witness=s.witness_times_s;
assert(sum(vecnorm(s.delta_v_km_s,2,2))<=c.search_max_dv_km_s*c.rejected_proposal_factor, ...
 'ctocscreen:v3:searchCostRejected','User search DV threshold exceeded; not a physical infeasibility.');
applied=false; actual=node;
for stop=unique([departure finish])
 if start==departure&&~applied
  actual.state=xx.'; actual.t=departure; actual.visited=visited;
  xx(4:6)=xx(4:6)+dv; applied=true;
 end
 if stop>start
  if isempty(arc)
   [xx,~,sol]=ctocscreen.v3Arc(xx,start,stop,eph.model,c);
   h=ctocscreen.v3Height(sol,eph.model,c);
  else
   sol=arc.solution;
   assert(departure==node.t&&start==departure&&stop==finish ...
    &&sol.x(1)==start&&sol.x(end)==stop&&norm(sol.y(1:6,1)-xx)<1e-10, ...
    'ctocscreen:v3:arcMismatch','Cached arc must start at this actual post-impulse state.');
   xx=sol.y(1:6,end); h=arc.height;
  end
  assert(h.passed,'ctocscreen:v3:actionHeight','Action violates height.');
  hints=nan(35,1); if target>0, hints(target)=finish; end
  if ~isempty(arc)&&isfield(arc,'scan')
   scan=arc.scan;
  else
   scan=ctocscreen.v3ScanArc(sol,@(ids,t,mode)ctocscreen.v3QueryTargets(eph,ids,t,mode),c,hints,find(~visited));
  end
  added=scan.visited&~visited; witness(added)=scan.time_s(added); visited=visited|scan.visited;
 end
 start=stop;
end
s.witness_times_s=witness; s=ctocscreen.v3Normalize(s,eph.model);
child=node; child.schedule=s; child.state=xx.'; child.t=finish; child.visited=visited;
child.last_gain=sum(visited)-sum(node.visited); child.J=sum(vecnorm(s.delta_v_km_s,2,2));
child.policy_rejected=child.J>c.search_max_dv_km_s;
plane=ctocscreen.v3PlaneChange(actual.state,dv,c);
penalty=0; if isfield(node,'inclination_penalty'), penalty=node.inclination_penalty; end
child.inclination_penalty=penalty+plane.penalty;
child.last_action=struct('origin',struct('state',actual.state,'t',departure, ...
 'visited',actual.visited,'J',node.J),'dv',dv,'dt',finish-departure,'target',target,'plane',plane);
kind=1; if any(dv~=0), kind=2; end; if target>0, kind=3; end
child.keys{end+1}=ctocscreen.v3StateKey(actual,kind,target,finish-departure,c);
left=find(~visited); estimate=0;
if ~c.cost_guidance_enabled&&~isempty(left)&&finish<eph.model.horizon_s
 lookahead=min(max(c.action_min_duration_s,min(c.action_durations_s)),eph.model.horizon_s-finish);
 if isempty(arc)
  future=ctocscreen.v3Arc(xx,finish,finish+lookahead,eph.model,c);
  targets=ctocscreen.v3QueryTargets(eph,left,finish+lookahead);
  miss=min(vecnorm(targets-future(1:3).',2,2));
 else
  % Cheap ranking only: finalizing an existing arc starts no new propagation.
  targets=ctocscreen.v3QueryTargets(eph,left,finish);
  miss=min(vecnorm(targets-xx(1:3).',2,2));
 end
 estimate=miss/lookahead*sqrt(numel(left));
end
child.estimate=child.J+estimate;
if c.cost_guidance_enabled
 [radial,shape]=ctocscreen.v3ContinuationEstimate(xx,visited,eph);
 child.estimate=child.J+c.arrival_shape_weight*shape+c.radial_tour_weight*radial;
end
child.estimate=child.estimate+c.plane_penalty_km_s*child.inclination_penalty;
end
