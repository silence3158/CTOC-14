function [out,tr]=replay(q,eph,c,independent,parent)
%REPLAY Execute fixed Cartesian controls; independent mode never uses a prefix.
if nargin<4, independent=false; end
if nargin<5, parent=[]; end
m=eph.model;
out=struct('passed',false,'status','invalid','visit_count',0,'distance_km',inf(35,1), ...
 'witness_times_s',nan(35,1),'total_dv_km_s',Inf,'height_passed',false, ...
 'min_altitude_lower_km',Inf,'sampled_min_altitude_km',Inf,'initial_passed',false, ...
 'inclination_penalty',0,'inclination_changes_deg',[],'plane_angles_deg',[], ...
 'independent',independent,'official_alignment_verified',false, ...
 'failure_id','','failure_reason','','reused_prefix',false);
tr=struct('q',q,'times',[],'pre',zeros(6,0),'post',zeros(6,0), ...
 'arcs',{{}},'arc_start',[],'arc_end',[],'heights',{{}},'scans',{{}},'key','');
out.dataset_kind='competition_targets_nominal_model';
if isfield(eph,'synthetic')&&eph.synthetic, out.dataset_kind='synthetic_test_only'; end
try
 q=ctocscreen.v4.canonical(q,m); tr.q=q; tr.key=ctocscreen.v4.controlKey(q);
 o=ctocscreen.v4.initialOrbit(q.x0,m); out.initial_passed=o.passed; out.initial_orbit=o;
 tq=@(ids,t,mode)ctocscreen.v3QueryTargets(eph,ids,t,mode);
 if independent&&q.T>0
  opt=odeset('RelTol',c.verify_reltol,'AbsTol',repmat([1e-12*ones(3,1);1e-15*ones(3,1)],35,1),'MaxStep',60);
  ts=ode89(@targetRhs,[0 q.T],reshape(eph.states0.',[],1),opt);
  assert(ts.x(end)>=q.T,'ctocscreen:v4:targetPropagation','Target propagation incomplete.');
  tq=@targetQuery;
 end
 start=0; x=q.x0; height=true; offset=0;
 if ~independent&&~isempty(parent)&&canReuse(q,parent.q)&&parent.actual.height_passed
  start=parent.q.T; x=parent.trace.post(:,end); out=parent.actual;
  out.independent=false; out.reused_prefix=true; out.passed=false;
  tr=parent.trace; tr.q=q; tr.key=ctocscreen.v4.controlKey(q);
  offset=numel(tr.times)-1;
 elseif ~independent&&~isempty(parent)
  % Same controls through the parent's last impulse: keep the arcs before it,
  % recompute only from there. Identical knots give the from-epoch computation.
  j=truncatedReuse(q,parent);
  if j>0
   keep=1:j-1; tr=parent.trace; tr.q=q; tr.key=ctocscreen.v4.controlKey(q);
   tr.times=tr.times(1:j); tr.pre=tr.pre(:,1:j); tr.post=tr.post(:,1:j);
   tr.arcs=tr.arcs(keep); tr.arc_start=tr.arc_start(keep); tr.arc_end=tr.arc_end(keep);
   tr.heights=tr.heights(keep); tr.scans=tr.scans(keep);
   for k=keep
    h=tr.heights{k}; height=height&&h.passed;
    out.min_altitude_lower_km=min(out.min_altitude_lower_km,h.lower_bound_altitude_km);
    out.sampled_min_altitude_km=min(out.sampled_min_altitude_km,h.sampled_min_altitude_km);
    % Scans read witness hints; rescan an arc whose relevant hints changed.
    a0=tr.arc_start(k); a1=tr.arc_end(k); w1=q.witness(:); w0=parent.trace.q.witness(:);
    inside=(w1>=a0&w1<=a1)|(w0>=a0&w0<=a1);
    if ~isequaln(w1(inside),w0(inside))
     visits=ctocscreen.v3ScanArc(tr.arcs{k},tq,c,q.witness);
     tr.scans{k}=struct('distance_km',visits.distance_km,'time_s',visits.time_s);
    end
    better=tr.scans{k}.distance_km<out.distance_km;
    out.distance_km(better)=tr.scans{k}.distance_km(better);
    out.witness_times_s(better)=tr.scans{k}.time_s(better);
   end
   start=tr.times(j); x=tr.pre(:,j); offset=j-1; out.reused_prefix=true;
  end
 end
 knots=unique([start;q.tau(q.tau>=start);q.T]);
 for k=1:numel(knots)
  t=knots(k); b=find(q.tau==t); before=x;
  % A reused endpoint already contains any original endpoint pulse.
  if ~isempty(b), x(4:6)=x(4:6)+sum(q.u(b,:),1).'; end
  ix=offset+k; tr.times(ix)=t; tr.pre(:,ix)=before; tr.post(:,ix)=x;
  if k==numel(knots), break; end
  % Construction commits the accurate six-state flow also used at acceptance.
  % Independent mode still starts at the epoch and integrates all targets anew.
  [x,~,sol]=ctocscreen.v3Arc(x,t,knots(k+1),m,c,false,true);
  h=ctocscreen.v3Height(sol,m,c); height=height&&h.passed;
  tr.arcs{end+1}=sol; tr.arc_start(end+1)=t; tr.arc_end(end+1)=knots(k+1);
  tr.heights{end+1}=h;
  out.min_altitude_lower_km=min(out.min_altitude_lower_km,h.lower_bound_altitude_km);
  out.sampled_min_altitude_km=min(out.sampled_min_altitude_km,h.sampled_min_altitude_km);
  visits=ctocscreen.v3ScanArc(sol,tq,c,q.witness);
  tr.scans{end+1}=struct('distance_km',visits.distance_km,'time_s',visits.time_s);
  better=visits.distance_km<out.distance_km;
  out.distance_km(better)=visits.distance_km(better);
  out.witness_times_s(better)=visits.time_s(better);
 end
 if q.T==0
  rt=eph.states0(:,1:3); out.distance_km=vecnorm(rt-q.x0(1:3).',2,2);
  out.witness_times_s=zeros(35,1); out.min_altitude_lower_km=norm(q.x0(1:3))-m.re;
  out.sampled_min_altitude_km=out.min_altitude_lower_km; height=out.min_altitude_lower_km>=200;
 end
 out.total_dv_km_s=sum(vecnorm(q.u,2,2)); out.height_passed=height;
 out.visit_count=sum(out.distance_km<=1); out.final_state=x;
 out.initial_passed=o.passed; out.initial_orbit=o;
 [out.inclination_penalty,out.inclination_changes_deg,out.plane_angles_deg]=planePenalty(q,tr,c);
 out.passed=o.passed&&height&&out.visit_count==35;
 out.status='partial'; if out.passed, out.status='complete'; end
 if ~height, out.status='height_failed'; end
 if ~o.passed, out.status='initial_orbit_failed'; end
 out.validation_level='nominal_j2_screened';
 if independent, out.validation_level='nominal_j2_independent_checked'; end
 if independent&&out.passed, out.validation_level='nominal_j2_independent_verified'; end
 out.official_alignment_verified=m.official_alignment_verified;
 out.event_scope='positive witnesses for all 35 targets; no proof of absence';
catch err
 out.status='propagation_failure'; out.failure_id=err.identifier; out.failure_reason=err.message;
 out.passed=false;
end
 function dx=targetRhs(t,y)
  s=reshape(y,6,35).'; a=ctocscreen.v3ReferenceForce(t,s(:,1:3),m);
  dx=reshape([s(:,4:6),a].',[],1);
 end
 function [r,v]=targetQuery(ids,t,mode)
  if isempty(ids), ids=1:35; end
  ids=ids(:); t=t(:); ni=numel(ids); nt=numel(t);
  if strcmp(mode,'grid'), ids=repmat(ids,nt,1); t=repelem(t,ni);
  elseif isscalar(ids), ids=repmat(ids,nt,1);
  elseif isscalar(t), t=repmat(t,ni,1); end
  [ut,~,map]=unique(t); yy=deval(ts,ut); r=zeros(numel(ids),3); v=r;
  for j=1:numel(ids)
   rows=6*(ids(j)-1); r(j,:)=yy(rows+(1:3),map(j)).'; v(j,:)=yy(rows+(4:6),map(j)).';
  end
  if strcmp(mode,'grid')
   r=permute(reshape(r,ni,nt,3),[1 3 2]); v=permute(reshape(v,ni,nt,3),[1 3 2]);
  end
 end
end
function j=truncatedReuse(q,parent)
% Trace index of the latest parent impulse knot t* such that the child also
% has an impulse at t* and identical controls before it; 0 when none exists.
% The child's knots before t* then match, so reuse equals the epoch replay.
j=0; p=parent.q; tr=parent.trace;
if isempty(p.tau)||~isfield(tr,'scans')||numel(tr.scans)~=numel(tr.arcs), return; end
if ~isequal(q.x0,p.x0)||~strcmp(tr.key,ctocscreen.v4.controlKey(p)), return; end
for k=numel(p.tau):-1:1
 t=p.tau(k); before=q.tau<t;
 if ~any(q.tau==t)||sum(before)~=k-1, continue; end
 if ~isequal(q.tau(before),p.tau(1:k-1))||~isequal(q.u(before,:),p.u(1:k-1,:)), continue; end
 i=find(tr.times==t,1);
 if ~isempty(i)&&numel(tr.arc_end)>=i-1&&all(tr.arc_end(1:i-1)<=t), j=i; end
 return
end
end
function yes=canReuse(q,p)
n=numel(p.tau);
yes=q.T>=p.T&&p.T>0&&isequal(q.x0,p.x0)&&numel(q.tau)>=n ...
 &&isequal(q.tau(1:n),p.tau)&&isequal(q.u(1:n,:),p.u) ...
 &&all(q.tau(n+1:end)>=p.T)&&all(p.tau<p.T) ...
 &&any(q.tau(n+1:end)==p.T);
end
function [pen,inc,plane]=planePenalty(q,tr,c)
inc=zeros(numel(q.tau),1); plane=inc;
for k=1:numel(q.tau)
 j=find(tr.times==q.tau(k),1); r=tr.pre(1:3,j); v=tr.pre(4:6,j);
 h=cross(r,v); g=cross(r,v+q.u(k,:).');
 assert(norm(h)>0&&norm(g)>0,'ctocscreen:v4:plane','Undefined orbit plane.');
 h=h/norm(h); g=g/norm(g);
 inc(k)=abs(acosd(max(-1,min(1,h(3))))-acosd(max(-1,min(1,g(3)))));
 plane(k)=acosd(max(-1,min(1,dot(h,g))));
end
pen=sum((max(0,inc-c.plane_threshold_deg)/c.plane_threshold_deg).^2);
end
