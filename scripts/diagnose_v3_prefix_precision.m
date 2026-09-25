function report=diagnose_v3_prefix_precision()
%DIAGNOSE_V3_PREFIX_PRECISION Read an old failure only to measure integration error.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
d=load(fullfile(sim,'runs/v3/search/v3_cold_1440_seed888_01/checkpoint.mat'),'state');
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
s=d.state.best_partial.schedule; c=d.state.config;
clock=tic; reference=integrate(s,eph.model,c,true,false);
report=struct('purpose','diagnostic_only_not_a_search_seed','reference_seconds',toc(clock),'cases',{{}});
for tolerance=[1e-10 1e-12 3e-14]
 for variational=[false true]
  local=c; local.shooting_reltol=tolerance; clock=tic;
  trace=integrate(s,eph.model,local,false,variational); elapsed=toc(clock); dr=0; dv=0;
  for k=1:numel(trace)
   a=trace{k}; times=linspace(a.x(1),a.x(end),max(3,ceil(diff(a.x([1 end]))/300)+1));
   x=deval(a,times); y=deval(reference{k},times);
   dr=max(dr,max(vecnorm(x(1:3,:)-y(1:3,:),2,1)));
   dv=max(dv,max(vecnorm(x(4:6,:)-y(4:6,:),2,1)));
  end
  item=struct('reltol',tolerance,'variational',variational,'seconds',elapsed,'position_difference_km',dr,'velocity_difference_km_s',dv);
  report.cases{end+1}=item;
  fprintf('rt %.3g stm %d wall %.3f dr %.9g km dv %.9g km/s\n',tolerance,variational,elapsed,dr,dv);
 end
end
save(fullfile(sim,'runs/v3/development/prefix_precision_diagnosis.mat'),'report');
end

function arcs=integrate(s,m,c,independent,variational)
x=ctocscreen.initialState(s.initial_q,m.mu,m.re).'; times=unique([0;s.maneuver_times_s;s.duration_s]); arcs={};
for k=1:numel(times)-1
 b=find(s.maneuver_times_s==times(k));
 if ~isempty(b), x(4:6)=x(4:6)+s.delta_v_km_s(b,:).'; end
 [x,~,arcs{end+1}]=ctocscreen.v3Arc(x,times(k),times(k+1),m,c,variational,independent);
end
end
