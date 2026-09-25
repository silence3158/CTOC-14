function [path,stats]=solve(times,positions,vStart,mu,re,stop)
%SOLVE Epoch-pair DP adapted from Zhang's MIT-licensed reference algorithm.
% Upstream: zhong-zh15/Multi_Flyby_Dynamic_Programming, 22126db5374f6eb465965a804d1a69706aefe1c0.
% Extension: retain every enumerated Lambert branch instead of preselecting
% a rendezvous-cheapest branch. The graph is a two-body proposal only.
n=numel(times)-1; edges=cell(n,1); costs=cell(n,1); prev=cell(n,1);
stats=struct('pairs',0,'branches',0,'height_rejected',0,'elapsed_s',0);
timer=tic; policy=struct('max_revolutions',Inf,'endpoint_tol_km',.001);
problem=struct('mu_km3_s2',mu,'re_km',re);
for k=1:n
 E=zeros(0,8); % departure epoch index, arrival index, departure/arrival velocity
 for a=1:numel(times{k})
  for b=1:numel(times{k+1})
   if stop(), error('zhangdp:budget','Graph construction budget exhausted.'); end
   dt=times{k+1}(b)-times{k}(a);
   if dt<=1, continue; end % Recorded numerical search approximation, not task constraint.
   policy.max_revolutions=ctocscreen.v3LambertRevolutions(positions{k}(a,:),positions{k+1}(b,:),dt,mu,Inf);
   bs=ctocscreen.v3LambertBranches(positions{k}(a,:),positions{k+1}(b,:),dt,mu,policy);
   stats.pairs=stats.pairs+1;
   for j=1:numel(bs)
    h=ctocscreen.checkArc([positions{k}(a,:),bs(j).v_depart],dt,problem);
    if ~strcmp(h.status,'ok')||h.min_altitude_km<200
     stats.height_rejected=stats.height_rejected+1; continue;
    end
    E(end+1,:)=[a,b,bs(j).v_depart,bs(j).v_arrive]; %#ok<AGROW>
   end
  end
 end
 edges{k}=E; stats.branches=stats.branches+size(E,1);
 if isempty(E), error('zhangdp:noPath','No admissible proposal arcs at stage %d.',k); end
 if k==1
  costs{k}=vecnorm(E(:,3:5)-vStart,2,2); prev{k}=zeros(size(E,1),1);
 else
  costs{k}=inf(size(E,1),1); prev{k}=zeros(size(E,1),1);
  last=edges{k-1}; lastCosts=costs{k-1};
  for a=1:numel(times{k})
   incoming=find(last(:,2)==a & isfinite(lastCosts)); outgoing=find(E(:,1)==a);
   if isempty(incoming), continue; end
   for j=outgoing.'
    [value,ix]=min(lastCosts(incoming)+vecnorm(last(incoming,6:8)-E(j,3:5),2,2));
    costs{k}(j)=value; prev{k}(j)=incoming(ix);
   end
  end
 end
 if ~any(isfinite(costs{k})), error('zhangdp:noPath','No connected path at stage %d.',k); end
end
[value,index]=min(costs{n}); selected=zeros(n,8);
for k=n:-1:1
 selected(k,:)=edges{k}(index,:); index=prev{k}(index);
end
arrival=zeros(n,1); departure=arrival; r0=zeros(n,3); r1=r0;
for k=1:n
 departure(k)=times{k}(selected(k,1)); arrival(k)=times{k+1}(selected(k,2));
 r0(k,:)=positions{k}(selected(k,1),:); r1(k,:)=positions{k+1}(selected(k,2),:);
end
path=struct('nominal_dv_km_s',value,'departure_times_s',departure, ...
 'arrival_times_s',arrival,'v_depart',selected(:,3:5),'v_arrive',selected(:,6:8), ...
 'r_depart',r0,'r_arrive',r1,'kind','fixed_order_visit_bound_two_body_proposal');
stats.elapsed_s=toc(timer);
end
