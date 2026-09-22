function f=replayMultiFragment(xin,tin,ids,counts,z,p)
%REPLAYMULTIFRAGMENT Generic burn/visit topology, exact natural propagation.
% counts=[0 3]: two burns, then three visits with no intervening impulse.
tokens=[]; id=0;
for m=1:numel(counts)
 tokens(end+1)=0; %#ok<AGROW>
 for j=1:counts(m), id=id+1; tokens(end+1)=id; end %#ok<AGROW>
end
L=numel(tokens); M=numel(counts); dt=z(1:L); u=reshape(z(L+1:end),3,M)';
f=struct('passed',false,'status','invalid','counts',counts,'state_in',xin,'t_in_s',tin, ...
 'maneuver_times_s',zeros(M,1),'delta_v_km_s',u,'event_target_ids',ids(:), ...
 'event_times_s',zeros(numel(ids),1),'event_states',zeros(numel(ids),6), ...
 'position_residual_km',zeros(numel(ids),3),'min_altitude_km',Inf);
assert(sum(counts)==numel(ids)&&numel(z)==L+3*M&&all(dt>=0));
assert(all(dt(2:end)>0)&&counts(end)>0&&all(counts>=0)&all(counts==floor(counts)));
x=xin; t=tin; m=0;
for j=1:L
 a=ctocscreen.checkArc(x,dt(j),p); assert(strcmp(a.status,'ok'));
 f.min_altitude_km=min(f.min_altitude_km,a.min_altitude_km);
 [x,info]=ctocscreen.propagateTwoBody(x,dt(j),p.mu_km3_s2); assert(strcmp(info.status,'ok'));
 t=t+dt(j);
 if tokens(j)==0
  m=m+1; f.maneuver_times_s(m)=t; x(4:6)=x(4:6)+u(m,:);
 else
  e=tokens(j); target=ctocscreen.targetStates(p,ids(e),t);
  f.event_times_s(e)=t; f.event_states(e,:)=x;
  f.position_residual_km(e,:)=x(1:3)-target(1:3);
 end
end
f.event_distances_km=vecnorm(f.position_residual_km,2,2);
f.total_dv_km_s=sum(vecnorm(u,2,2)); f.state_out=x; f.t_out_s=t;
f.passed=all(f.event_distances_km<=.05)&&f.min_altitude_km>=200&&t<=p.horizon_s;
f.status='fragment_infeasible'; if f.passed, f.status='fragment_screened'; end
end
