function [children,failures,report]=v3GuidedChildren(node,actions,eph,c)
%V3GUIDEDCHILDREN Finish successful arcs before diversity pruning, even at deadline.
children={}; failures={};
report=struct('processed',0,'generated',0,'retained',0,'reused_arcs',0);
for k=1:numel(actions)
 a=actions{k}; report.processed=report.processed+1;
 try
  children{end+1}=ctocscreen.v3ApplyAction(node,node.t,a.dv,node.t+a.dt,a.target,eph,c,a.arc);
  report.reused_arcs=report.reused_arcs+1;
 catch err
  failures{end+1}=struct('id',err.identifier,'message',err.message);
 end
end
report.generated=numel(children);
children=ctocscreen.v3SelectBeam(children,c.guided_keep,c);
report.retained=numel(children);
end
