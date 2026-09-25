function report=v3PolicyFeedback(table,node,target,dt,dv,c,rejectedCost,rejectedComplete)
%V3POLICYFEEDBACK Bounded negative evidence from actual impulses, not ODE failure.
% Graded credit assignment (BWAS "only the responsible ant is updated", adapted to
% a continuous spend): inside-budget burns are never charged, the crossing burn
% carries the strongest share, later burns of an overspent candidate carry their
% own share. Coasting cannot refund DV already spent.
plane=ctocscreen.v3PlaneChange(node.state,dv,c);
spend=node.J+norm(dv);
total=spend; excess=max(0,total-c.search_max_dv_km_s);
covered=node.visited; if target>0, covered(target)=true; end
complete=all(covered);
explicitRejection=nargin>=7&&isfinite(rejectedCost);
if explicitRejection
 complete=true; if nargin>=8, complete=rejectedComplete; end
 total=rejectedCost; excess=max(0,rejectedCost-c.search_max_dv_km_s);
 % Consistency guard: rejectedCost must be the FULL candidate cost including
 % this burn. Passing a prefix cost would silently invert the severity ranking.
 assert(rejectedCost>=spend-1e-9,'ctocscreen:v3:policyFeedback', ...
  'rejectedCost must include the current burn (prefix cost inverted the ranking).');
elseif ~(complete&&isfinite(c.search_max_dv_km_s))
 excess=0;
end
% Graded credit assignment (BWAS "only the responsible ant is updated", adapted
% to a continuous spend): a burn that is still inside the budget is never charged;
% the burn that crosses the cap carries the strongest share; later burns of an
% overspent candidate carry their own share so mid-course waste stays learnable.
% Coasting cannot refund DV already spent.
if isfinite(c.search_max_dv_km_s)
 burned=node.J;
 if burned+norm(dv)<=c.search_max_dv_km_s
  credit=0;
 elseif burned<c.search_max_dv_km_s
  credit=1;
 else
  credit=norm(dv)/max(spend,eps);
 end
else
 credit=0;
end
costPenalty=0;
if excess>0&&credit>0
 % Monotone in the overspend and scaled by the responsible share of the real
 % cost, so the stored penalty grows with severity. The global cap is only the
 % unit here; it is not the accepted threshold and does not clip the ranking.
 share=norm(dv)/max(total,eps);
 costPenalty=excess/max(c.search_max_dv_km_s,eps)*share*credit;
end
% The 1/(1+penalty) query then makes a larger penalty a strictly smaller weight,
% which is what actually steers the search away from the expensive structure.
penalty=min(c.policy_penalty_ceiling,(plane.penalty+costPenalty)*c.policy_feedback_gain);
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
burned=spend-norm(dv);
crossing=credit>0&&burned<c.search_max_dv_km_s;
report=struct('total_dv_km_s',total,'complete_candidate',complete, ...
 'rejected_cost',(complete||explicitRejection)&&excess>0, ...
 'rejected_prefix',explicitRejection&&~complete&&excess>0, ...
 'crossing_burn',crossing,'credit',credit,'cost_penalty',costPenalty, ...
 'inclination_change_deg',plane.inclination_change_deg,'plane_rotation_deg',plane.plane_rotation_deg, ...
 'penalty',penalty,'new_penalty',increment,'policy_key',key);
end
