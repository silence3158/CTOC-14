function p=v3ShootingProblem(s,eph,c,cached)
%V3SHOOTINGPROBLEM Scaled sparse multiple shooting with free event times.
% Event ordering is fixed for one NLP only; the outer mutation changes it.
s=ctocscreen.v3Normalize(s,eph.model); m=eph.model;
ids=c.joint_target_ids(:); K=numel(ids); w=s.witness_times_s(ids);
if isfield(s,'visit_plan_times_s'), w=s.visit_plan_times_s(ids); end
missing=isnan(w);
if any(missing)
 nearest=ctocscreen.v3Replay(s,eph,c,false);
 assert(~strcmp(nearest.status,'propagation_failure'),'ctocscreen:v3:seedReplay','Cannot initialize from failed propagation.');
 w(missing)=nearest.witness_times_s(ids(missing));
end
base=unique([0;s.maneuver_times_s;w;s.duration_s]); extra=[];
for segment=1:numel(base)-1
 tt=linspace(base(segment),base(segment+1),max(2,ceil(diff(base(segment:segment+1))/c.max_shooting_arc_s)+1));
 extra=[extra;reshape(tt(2:end-1),[],1)]; %#ok<AGROW>
end
% Keep coincident visits/burns as separate movable events, including tau=0.
% Only the first epoch is fixed. Zero-time arcs carry exact identity defects.
M=numel(s.maneuver_times_s);
[t,eventOrder]=sort([0;s.maneuver_times_s;w;s.duration_s;extra]); N=numel(t);
[~,bn]=ismember(1+(1:M),eventOrder); [~,wn]=ismember(1+M+(1:K),eventOrder);
sx=[1e4;1e4;1e4;1;1;1]; sq=[10;.001;.001;1;1;1]; qbase=[m.re+600;0;0;0;0;0]; st=86400;
iq=1:6; it=[0,6+(1:N-1)]; ix=reshape(6+N-1+(1:6*N),6,N);
idv=reshape(max(ix(:))+(1:3*M),3,M); is=max(ix(:))+3*M+(1:M);
n=6+N-1+6*N+4*M; z=zeros(n,1);
z(iq)=(s.initial_q(:)-qbase)./sq; z(it(2:end))=t(2:end)/st;
if nargin>=4&&~isempty(cached)
 assert(strcmp(cached.schedule_key,ctocscreen.v3ScheduleKey(s)), ...
  'ctocscreen:v3:traceMismatch','Initialization trace must describe these exact controls.');
 replayed=cached.replay; trace=cached.trace;
else
 [replayed,trace]=ctocscreen.v3Replay(s,eph,c,false);
end
assert(~strcmp(replayed.status,'propagation_failure'),'ctocscreen:v3:seedReplay','Seed propagation failed.');
for node=1:N
 knot=find(trace.times==t(node),1);
 if isempty(knot)
  arc=find(trace.times<t(node),1,'last'); x=deval(trace.arcs{arc},t(node)); x=x(1:6);
 else
  x=trace.states(knot,:).';
  % Replay knots are post-burn; visits ordered before a coincident burn
  % need its pre-burn state. Event indices remain independently movable.
  futureBurn=find(s.maneuver_times_s==t(node) & bn(:)>node);
  if ~isempty(futureBurn), x(4:6)=x(4:6)-sum(s.delta_v_km_s(futureBurn,:),1).'; end
 end
 z(ix(:,node))=x./sx;
end
z(idv(:))=reshape(s.delta_v_km_s.',[],1); z(is)=vecnorm(s.delta_v_km_s,2,2);
% All initial nodes are actual propagated states: no artificial target resets.
lb=-inf(n,1); ub=inf(n,1); lb(iq)=[-1;-1;-1;0;-Inf;-Inf]; ub(iq)=[1;1;1;pi;Inf;Inf];
lb(it(2:end))=0; ub(it(2:end))=m.horizon_s/st;
lb(ix(1:3,:))=-c.node_position_bound/1e4; ub(ix(1:3,:))=c.node_position_bound/1e4;
lb(ix(4:6,:))=-c.node_velocity_bound; ub(ix(4:6,:))=c.node_velocity_bound;
lb(idv(:))=-c.pulse_component_bound; ub(idv(:))=c.pulse_component_bound;
lb(is)=0; ub(is)=sqrt(3)*c.pulse_component_bound;
% These are expandable search bounds, not physical limits. Preserve the
% input trajectory rather than silently clipping nodes or pulse components.
expand=[ix(:);idv(:);is(:)];
expanded=sum(z(expand)<lb(expand)|z(expand)>ub(expand));
lb(expand)=min(lb(expand),z(expand)-1e-9); ub(expand)=max(ub(expand),z(expand)+1e-9);
B=sparse(N-2,n);
for event=2:N-1, B(event-1,it(event))=1; B(event-1,it(event+1))=-1; end
lastZ=[]; lastC=[]; lastE=[]; lastG=[]; lastGE=[];
p=struct('z0',min(ub,max(lb,z)),'lb',lb,'ub',ub,'A',B,'b',zeros(N-2,1), ...
 'objective',@objective,'constraints',@constraints,'decode',@decode,'connect',@connect, ...
 'encounter_linearization',@encounterLinearization, ...
 'node_count',N,'maneuver_count',M,'variable_count',n,'witness_nodes',wn, ...
 'time_indices',it,'state_indices',ix,'q_indices',iq,'active_target_ids',ids, ...
 'expanded_bound_count',expanded);
 function [f,g]=objective(v)
  f=sum(v(is)); g=zeros(n,1); g(is)=1;
 end
 function cand=decode(v)
  tt=[0;v(it(2:end))*st]; cand=s;
  cand.initial_q=(qbase+sq.*v(iq)).'; cand.maneuver_times_s=tt(bn);
  cand.delta_v_km_s=reshape(v(idv(:)),3,M).'; cand.duration_s=tt(end);
   cand.witness_times_s(ids)=tt(wn); cand.validation_level='proposal';
   if isfield(cand,'visit_plan_times_s'), cand.visit_plan_times_s(ids)=tt(wn); end
 end
 function [cand,report]=connect(v,budget,mode)
  if nargin<3, mode='free_events'; end
  % Retract solver event coordinates, never clip an input physical schedule.
  rawTimes=[0;v(it(2:end))*st];
  times=cummax(max(0,min(m.horizon_s,rawTimes)));
  v(it(2:end))=times(2:end)/st;
  cand=decode(v); X=reshape(v(ix(:)),6,N).*sx;
  cand.maneuver_times_s=times(bn); cand.duration_s=times(end);
  cand.witness_times_s(ids)=times(wn);
  if isfield(cand,'visit_plan_times_s'), cand.visit_plan_times_s(ids)=times(wn); end
  anchors=X(1:3,[bn(2:end),N]).';
  if M==0, anchors=zeros(0,3); end
  aligned=0;
  if strcmp(mode,'encounter_aligned')
   endpoints=[bn(2:end),N];
   for j=1:M
    match=find(abs(t(wn)-t(endpoints(j)))<=1e-3);
    if numel(match)~=1, continue; end
    id=ids(match); epoch=times(endpoints(j));
    anchors(j,:)=ctocscreen.v3QueryTargets(eph,id,epoch);
    cand.witness_times_s(id)=epoch;
    if isfield(cand,'visit_plan_times_s'), cand.visit_plan_times_s(id)=epoch; end
    aligned=aligned+1;
   end
  end
  merged=0;
  if M>1
   [burns,~,group]=unique(cand.maneuver_times_s); pulses=zeros(numel(burns),3);
   merged=M-numel(burns); ends=zeros(numel(burns),1);
   for j=1:M
    pulses(group(j),:)=pulses(group(j),:)+cand.delta_v_km_s(j,:);
    ends(group(j))=j;
   end
   cand.maneuver_times_s=burns; cand.delta_v_km_s=pulses; anchors=anchors(ends,:);
  end
  [cand,report]=ctocscreen.v3ConnectControls(cand,anchors,m,c,budget);
  report.time_projection_max_s=max(abs(times-rawTimes));
  report.raw_time_range_s=[min(rawTimes),max(rawTimes)];
  report.raw_min_event_gap_s=min(diff(rawTimes)); report.merged_burns=merged;
  report.projection_mode=mode; report.aligned_visits=aligned;
 end
 function [Avisit,bvisit]=encounterLinearization(v)
  % Inner box of the visit ball is a QP approximation, not a formal limit.
  % Unlike squared distance at the centre, vector position has full slope.
  tt=[0;v(it(2:end))*st]; X=reshape(v(ix(:)),6,N).*sx;
  J=sparse(3*K,n); residual=zeros(3*K,1);
  for j=1:K
   k=wn(j); rows=3*(j-1)+(1:3);
   [rt,vt]=ctocscreen.v3QueryTargets(eph,ids(j),tt(k));
   residual(rows)=(X(1:3,k)-rt.')/1e4;
   J(rows,ix(1:3,k))=eye(3);
   if k>1, J(rows,it(k))=-vt.'*st/1e4; end
  end
  halfwidth=c.search_radius_km/sqrt(3)/1e4;
  Avisit=[J;-J]; bvisit=[halfwidth-residual;halfwidth+residual];
 end
 function [C,E,G,GE]=constraints(v)
  if isequal(v,lastZ), C=lastC; E=lastE; G=lastG; GE=lastGE; return; end
  % Solvers may evaluate just outside time bounds. Keep dynamics finite;
  % bounds/ordering still reject those trial coordinates.
  tt=max(0,min(m.horizon_s,[0;v(it(2:end))*st])); X=reshape(v(ix(:)),6,N).*sx;
  D=reshape(v(idv(:)),3,M); q=qbase+sq.*v(iq);
  C=zeros(1+M+K+N+N-1,1); E=zeros(6*N,1);
  G=sparse(n,numel(C)); GE=sparse(n,numel(E));
  C(1)=v(2)^2+v(3)^2-(1-1e-8)^2; G(2:3,1)=2*v(2:3);
  for b=1:M
   col=1+b; C(col)=sum(D(:,b).^2)-v(is(b))^2;
   G(idv(:,b),col)=2*D(:,b); G(is(b),col)=-2*v(is(b));
  end
  y0=ctocscreen.initialState(q,m.mu,m.re).'; Jq=zeros(6,6);
  for j=1:6
   hq=1e-5; qp=q; qm=q; qp(j)=qp(j)+hq*sq(j); qm(j)=qm(j)-hq*sq(j);
   Jq(:,j)=(ctocscreen.initialState(qp,m.mu,m.re)-ctocscreen.initialState(qm,m.mu,m.re)).'/(2*hq);
  end
  defect=X(:,1)-y0; b=find(bn==1);
  if ~isempty(b), defect(4:6)=defect(4:6)-D(:,b); GE(idv(:,b),4:6)=-eye(3); end
  E(1:6)=defect./sx; GE(ix(:,1),1:6)=eye(6); GE(iq,1:6)=(-Jq./sx).';
  for k=1:N-1
   % Infeasible trial event ordering is clamped only for finite evaluation;
   % linear inequalities still reject it. Feasible iterates use exact times.
   ta=tt(k); tb=max(ta,tt(k+1));
   [yf,P,sol]=ctocscreen.v3Arc(X(:,k),ta,tb,m,c,true);
   b=find(bn==k+1); d=zeros(6,1); if ~isempty(b), d(4:6)=D(:,b); end
   cols=6*k+(1:6); E(cols)=(X(:,k+1)-yf-d)./sx;
   GE(ix(:,k),cols)=(-P.*sx.'./sx).'; GE(ix(:,k+1),cols)=eye(6);
   if ~isempty(b), GE(idv(:,b),cols(4:6))=-eye(3); end
   fa=[X(4:6,k);ctocscreen.v3Force(ta,X(1:3,k),m)];
   fb=[yf(4:6);ctocscreen.v3Force(tb,yf(1:3),m)];
   if k>1, GE(it(k),cols)=(P*fa./sx*st).'; end
   GE(it(k+1),cols)=(-fb./sx*st).';
   col=1+M+K+N+k; mid=(ta+tb)/2;
   if isempty(sol), ym=X(:,k); Pm=eye(6);
   else, zm=deval(sol,mid); ym=zm(1:6); Pm=reshape(zm(7:end),6,6); end
   nr=norm(ym(1:3)); grad=[-ym(1:3).'/nr,zeros(1,3)]/1e4;
   C(col)=(m.re+200+c.height_margin_km-nr)/1e4;
   G(ix(:,k),col)=(grad*Pm.*sx.').';
   fm=[ym(4:6);ctocscreen.v3Force(mid,ym(1:3),m)];
   if k>1, G(it(k),col)=grad*(0.5*fm-Pm*fa)*st; end
   G(it(k+1),col)=grad*(0.5*fm)*st;
  end
  for j=1:K
   k=wn(j); [rt,vt]=ctocscreen.v3QueryTargets(eph,ids(j),tt(k)); dr=X(1:3,k)-rt.';
   % Distance and dynamics defects share the same position scale. Squared
   % misses otherwise dominate restoration and pull auxiliary nodes apart.
   radius=c.search_radius_km; distance=sqrt(dot(dr,dr)+radius^2); direction=dr/distance;
   col=1+M+j; C(col)=(distance-sqrt(2)*radius)/1e4;
   G(ix(1:3,k),col)=direction.*sx(1:3)/1e4;
   if k>1, G(it(k),col)=-dot(direction,vt)*st/1e4; end
  end
  for k=1:N
   col=1+M+K+k; nr=norm(X(1:3,k)); C(col)=(m.re+200+c.height_margin_km-nr)/1e4;
   G(ix(1:3,k),col)=-X(1:3,k)/nr;
  end
  lastZ=v; lastC=C; lastE=E; lastG=G; lastGE=GE;
 end
end
