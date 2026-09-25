function [candidate,name,report]=v3ConnectedMutation(original,eph,c,stream,budget,operator)
%V3CONNECTEDMUTATION Construct changed topology using all planned arc visits.
clock=tic;
if nargin<6
 [candidate,name]=ctocscreen.v3Mutate(original,eph,c,stream);
else
 [candidate,name]=ctocscreen.v3Mutate(original,eph,c,stream,operator);
end
candidate.visit_plan_times_s=candidate.witness_times_s;
[candidate,report]=ctocscreen.v3FitVisitArcs(candidate,eph,c,budget);
report.difference=ctocscreen.v3MutationDifference(original,candidate,c,eph.model);
report.constructed=report.passed;
report.visit_assignment_changed=~isequaln(candidate.visit_plan_times_s,original.witness_times_s);
report.passed=report.passed&&(report.difference.novel||report.visit_assignment_changed);
if report.constructed&&~report.passed
 report.reason='ctocscreen:v3:unchangedMutation: Correction returned to the parent trajectory.';
end
report.full_feasibility_verified=false;
name=['connected_' name]; report.elapsed_s=toc(clock);
end
