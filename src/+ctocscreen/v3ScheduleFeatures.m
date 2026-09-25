function f=v3ScheduleFeatures(s,eph,c,trace)
%V3SCHEDULEFEATURES Numerical diversity descriptor, never a feasibility test.
if nargin<4
 [r,trace]=ctocscreen.v3Replay(s,eph,c,false);
 assert(~strcmp(r.status,'propagation_failure'),'ctocscreen:v3:features','Cannot describe failed propagation.');
end
f=struct('physical_key',ctocscreen.v3ScheduleKey(s),'duration_s',s.duration_s, ...
 'initial_state',ctocscreen.initialState(s.initial_q,eph.model.mu,eph.model.re), ...
 'witness_times_s',s.witness_times_s,'states',zeros(c.diversity_samples,6));
keep=vecnorm(s.delta_v_km_s,2,2)>c.diversity_pulse_floor_km_s;
f.maneuver_times_s=s.maneuver_times_s(keep); f.delta_v_km_s=s.delta_v_km_s(keep,:);
[~,order]=sort(s.witness_times_s);
root=sprintf('%.17g/',s.initial_q); if isfield(s,'root_key'), root=s.root_key; end
f.family=[root '|order=' sprintf('%d/',order) '|burns=' num2str(sum(keep))];
tt=linspace(0,s.duration_s,c.diversity_samples);
for k=1:numel(tt)
 knot=find(trace.times==tt(k),1);
 if isempty(knot)
  arc=find(trace.times<tt(k),1,'last'); x=deval(trace.arcs{arc},tt(k)); f.states(k,:)=x(1:6).';
 else
  f.states(k,:)=trace.states(knot,:);
 end
end
end
