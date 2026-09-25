function [radial,shape]=v3ContinuationEstimate(state,visited,eph,ranges)
%V3CONTINUATIONESTIMATE Apsidal expansion surrogate, never a J2 lower bound.
% Uses actual velocity; a flyby need not circularize or match target velocity.
radial=0; shape=0; ids=find(~visited);
if isempty(ids), return; end
mu=eph.model.mu;
if nargin<4
 [peri,apo]=ctocscreen.v3OsculatingApses(eph.states0(ids,:),mu);
 ranges=[min(apo),max(peri)];
end
requiredPeri=ranges(1); requiredApo=ranges(2);
state=state(:).'; r=norm(state(1:3)); [rp,ra]=ctocscreen.v3OsculatingApses(state,mu);
if ~all(isfinite([rp ra requiredPeri requiredApo]))||rp<=0||ra<=0, return; end
p=min(rp,requiredPeri); a=max(ra,requiredApo);
% Compare the two orders of changing apsides. Waiting/phase are omitted.
outwardFirst=abs(vperi(rp,a)-vperi(rp,ra))+abs(vapo(p,a)-vapo(rp,a));
inwardFirst=abs(vapo(p,ra)-vapo(rp,ra))+abs(vperi(p,a)-vperi(p,ra));
radial=min(outwardFirst,inwardFirst);
if requiredApo<=1.1*requiredPeri
 % Circular coasting is useful only when remaining radial ranges overlap.
 vr=dot(state(1:3),state(4:6))/r;
 vt=norm(cross(state(1:3),state(4:6)))/r;
 shape=hypot(vr,vt-sqrt(mu/r));
end
 function v=vperi(rp,ra), v=sqrt(2*mu*ra/(rp*(rp+ra))); end
 function v=vapo(rp,ra), v=sqrt(2*mu*rp/(ra*(rp+ra))); end
end
