function events=scanEncounters(state,t0,duration,p,ids,step)
%SCANENCOUNTERS Candidate-only grid/local-minimum screen, not exhaustive proof.
if nargin<6, step=300; end
grid=linspace(0,duration,max(3,ceil(duration/step)+1));
d=zeros(numel(grid),numel(ids));
for k=1:numel(grid)
 x=ctocscreen.propagateTwoBody(state,grid(k),p.mu_km3_s2);
 y=ctocscreen.targetStates(p,ids,t0+grid(k));
 d(k,:)=vecnorm(y(:,1:3)-x(1:3),2,2)';
end
events=zeros(numel(ids),3);
for n=1:numel(ids)
 [best,j]=min(d(:,n)); when=grid(j);
 loc=find(d(2:end-1,n)<=d(1:end-2,n)&d(2:end-1,n)<=d(3:end,n))+1;
 for j=loc'
  [tt,dd]=fminbnd(@distance,grid(j-1),grid(j+1),optimset('Display','off','TolX',1e-5));
  if dd<best, best=dd; when=tt; end
 end
 events(n,:)=[ids(n) t0+when best];
end
events=sortrows(events,3);
 function dd=distance(dt)
  x=ctocscreen.propagateTwoBody(state,dt,p.mu_km3_s2);
  y=ctocscreen.targetStates(p,ids(n),t0+dt); dd=norm(x(1:3)-y(1:3));
 end
end
