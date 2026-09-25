function [keys,actions]=v3TraceExperience(s,trace,c)
%V3TRACEEXPERIENCE Credit actual pre-burn states from an existing matching replay.
keys={}; actions={};
for k=1:numel(trace.times)
 ta=trace.times(k); if ta>s.duration_s, break; end
 tb=min(s.duration_s,trace.times(min(k+1,numel(trace.times)))); x=trace.states(k,:);
 burn=find(s.maneuver_times_s==ta,1); kind=1;
 if ta==s.duration_s&&isempty(burn), continue; end
 if ~isempty(burn), x(4:6)=x(4:6)-s.delta_v_km_s(burn,:); kind=2; end
 node=struct('state',x,'t',ta,'visited',s.witness_times_s<=ta);
 ids=find(s.witness_times_s>ta & s.witness_times_s<=tb);
 keys{end+1}=ctocscreen.v3StateKey(node,kind,0,tb-ta,c);
 keys{end+1}=ctocscreen.v3PolicyKey(node,kind,0,tb-ta);
 if kind==2
  node.J=sum(vecnorm(s.delta_v_km_s(1:burn-1,:),2,2));
  target=0; dt=tb-ta;
  if ~isempty(ids), [~,first]=min(s.witness_times_s(ids)); target=ids(first); dt=s.witness_times_s(target)-ta; end
  actions{end+1}=struct('origin',node,'target',target,'dt',dt,'dv',s.delta_v_km_s(burn,:));
  for target=ids(:).'
   keys{end+1}=ctocscreen.v3StateKey(node,3,target,s.witness_times_s(target)-ta,c);
   keys{end+1}=ctocscreen.v3PolicyKey(node,3,target,s.witness_times_s(target)-ta);
  end
 end
end
keys=unique(keys);
end
