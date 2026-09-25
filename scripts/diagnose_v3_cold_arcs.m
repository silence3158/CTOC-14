function report=diagnose_v3_cold_arcs(label)
%DIAGNOSE_V3_COLD_ARCS Read-only fixed-control numerical sensitivity diagnosis.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
d=load(fullfile(sim,'runs/v3/search',label,'checkpoint.mat'),'state');
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
s=d.state.best_partial.schedule; c=d.state.config;
[a,ta]=ctocscreen.v3Replay(s,eph,c,false,'prefix_witnesses');
[b,tb]=ctocscreen.v3Replay(s,eph,c,true,'prefix_witnesses');
fine=c; fine.max_step_s=30;
[f,tf]=ctocscreen.v3Replay(s,eph,fine,true,'prefix_witnesses');
report=struct('purpose','diagnostic_only_not_initialization','screen',a,'reference',b,'fine',f,'arcs',[]);
for k=1:numel(ta.arcs)
 sa=ta.arcs{k}; sb=tb.arcs{k}; t0=sb.x(1); t1=sb.x(end);
 [local,P]=ctocscreen.v3Arc(sb.y(1:6,1),t0,t1,eph.model,c,true);
 scalar=ctocscreen.v3Arc(sb.y(1:6,1),t0,t1,eph.model,c);
 tt=linspace(t0,t1,100); xa=deval(sa,tt); xb=deval(sb,tt); xf=deval(tf.arcs{k},tt);
 row=[k,t1,max(vecnorm(xa(1:3,:)-xb(1:3,:),2,1)), ...
  max(vecnorm(xf(1:3,:)-xb(1:3,:),2,1)),norm(scalar(1:3)-sb.y(1:3,end)), ...
  norm(local(1:3)-sb.y(1:3,end)),norm(P),tb.heights{k}.sampled_min_altitude_km];
 report.arcs(end+1,:)=row;
 fprintf('arc%2d t%9.1f global %.3g refFine %.3g local %.3g STM %.3g gain %.3g h %.1f\n',row);
end
save(fullfile(sim,'runs/v3/search',label,'arc_diagnosis.mat'),'report');
end
