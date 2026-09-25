function report=v3PolicyFeedback(table,node,target,dt,dv,c,rejectedCost,rejectedComplete)
%V3POLICYFEEDBACK Bounded negative evidence from actual impulses, not ODE failure.
plane=ctocscreen.v3PlaneChange(node.state,dv,c);
total=node.J+norm(dv); excess=max(0,total-c.search_max_dv_km_s);
penalty=plane.penalty;
covered=node.visited; if target>0, covered(target)=true; end
% Overspent prefixes need earlier controls rebuilt; coasting cannot refund DV.
complete=all(covered);
explicitRejection=nargin>=7&&isfinite(rejectedCost);
if explicitRejection
 complete=true; if nargin>=8, complete=rejectedComplete; end
 excess=max(0,rejectedCost-c.search_max_dv_km_s);
 penalty=penalty+excess/c.search_max_dv_km_s*norm(dv)/max(rejectedCost,eps);
elseif complete&&isfinite(c.search_max_dv_km_s)
 penalty=penalty+excess/c.search_max_dv_km_s;
end
penalty=min(c.policy_penalty_ceiling,penalty*c.policy_feedback_gain);
kind=1; if any(dv~=0), kind=2; end; if target>0, kind=3; end
key=ctocscreen.v3PolicyKey(node,kind,target,dt);
identity=['observed|' ctocscreen.v3StateKey(node,kind,target,dt,c) '|' ...
 sprintf('%.17g/',[node.state(:);node.J;dt;dv(:)])];
old=0; if isKey(table,identity), old=table(identity); end
increment=max(0,penalty-old);
if increment>0
 table(identity)=penalty;
 ctocscreen.v3Pheromone('penalize',table,key,increment,c);
end
report=struct('total_dv_km_s',total,'complete_candidate',complete, ...
 'rejected_cost',(complete||explicitRejection)&&excess>0, ...
 'rejected_prefix',explicitRejection&&~complete&&excess>0, ...
 'inclination_change_deg',plane.inclination_change_deg,'plane_rotation_deg',plane.plane_rotation_deg, ...
 'penalty',penalty,'new_penalty',increment,'policy_key',key);
end
