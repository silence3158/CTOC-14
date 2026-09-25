function [path,stats]=delayed(times,order,eph,x,tFirst,waits,waitWidth,beam,stop,incumbent)
%DELAYED Finite-width extension of Zhang's state-dependent stage recurrence.
% Delayed departures require actual arrival velocity; epoch pairs cannot merge
% these labels. Keep beam labels per arrival epoch, as a search approximation.
% Coast and Lambert arcs are two-body proposals; final repair uses actual J2.
timer=tic; m=eph.model; n=numel(order); layers=cell(n,1);
if nargin<10, incumbent=[]; end
protectedParent=1;
labels=[tFirst,x(:).',0,0,0,zeros(1,3)];
% Columns: time, r(3), v(3), cost, parent index, departure time, departure v(3).
stats=struct('pairs',0,'branches',0,'coast_rejected',0,'height_rejected',0, ...
 'kept_per_stage',zeros(n,1),'elapsed_s',0,'beam_per_epoch',beam,'wait_width_s',waitWidth);
policy=struct('max_revolutions',Inf,'endpoint_tol_km',.001);
problem=struct('mu_km3_s2',m.mu,'re_km',m.re);
for k=1:n
 next=zeros(0,13); protectedRow=0; goals=ctocscreen.v3QueryTargets(eph,order(k),times{k},'pairs');
 if k==1, ws=0;
 elseif waitWidth==0, ws=waits(k);
 else, ws=unique(max(0,[0,waits(k)-waitWidth,waits(k),waits(k)+waitWidth])); end
 for p=1:size(labels,1)
  parent=labels(p,:);
  for w=ws
   if stop(), error('zhangdp:budget','Delayed construction budget exhausted at stage %d.',k); end
   td=parent(1)+w;
   if td>=max(times{k}), continue; end
   state=parent(2:7);
   if w>0
    h=ctocscreen.checkArc(state,w,problem);
    if ~strcmp(h.status,'ok')||h.min_altitude_km<200
     stats.coast_rejected=stats.coast_rejected+1; continue;
    end
    state=ctocscreen.propagateTwoBody(state,w,m.mu);
   end
   for q=1:numel(times{k})
    dt=times{k}(q)-td;
    if dt<=1, continue; end
    policy.max_revolutions=ctocscreen.v3LambertRevolutions(state(1:3),goals(q,:),dt,m.mu,Inf);
    bs=ctocscreen.v3LambertBranches(state(1:3),goals(q,:),dt,m.mu,policy);
    stats.pairs=stats.pairs+1;
    for b=1:numel(bs)
     h=ctocscreen.checkArc([state(1:3),bs(b).v_depart],dt,problem);
     if ~strcmp(h.status,'ok')||h.min_altitude_km<200
      stats.height_rejected=stats.height_rejected+1; continue;
     end
     cost=parent(8)+norm(bs(b).v_depart-state(4:6));
     next(end+1,:)=[times{k}(q),goals(q,:),bs(b).v_arrive,cost,p,td,bs(b).v_depart]; %#ok<AGROW>
     if ~isempty(incumbent)&&p==protectedParent&&abs(td-incumbent.departure_times_s(k))<1e-6 ...
       &&abs(times{k}(q)-incumbent.arrival_times_s(k))<1e-6 ...
       &&norm(bs(b).v_depart-incumbent.v_depart(k,:))<1e-6
      protectedRow=size(next,1);
     end
     stats.branches=stats.branches+1;
    end
   end
  end
 end
 if isempty(next), error('zhangdp:noPath','No delayed path at stage %d.',k); end
 keptIndices=zeros(0,1);
 for q=1:numel(times{k})
  ids=find(next(:,1)==times{k}(q)); [~,ix]=sort(next(ids,8)); selected=zeros(0,1);
  for j=reshape(ids(ix),1,[])
   % Remove only numerically identical endpoint velocities before beam cut.
   if isempty(selected)||all(vecnorm(next(selected,5:7)-next(j,5:7),2,2)>1e-9)
    selected(end+1,1)=j; %#ok<AGROW>
    if numel(selected)>=beam, break; end
   end
  end
  keptIndices=[keptIndices;selected]; %#ok<AGROW>
 end
 if ~isempty(incumbent)
  assert(protectedRow>0,'zhangdp:incumbentLost','Incumbent not reconstructed at stage %d.',k);
  if ~ismember(protectedRow,keptIndices), keptIndices(end+1)=protectedRow; end
  protectedParent=find(keptIndices==protectedRow,1);
 end
 labels=next(keptIndices,:); layers{k}=labels; stats.kept_per_stage(k)=size(labels,1);
 if mod(k,5)==0, fprintf('DELAYED stage %d prefix %.9f labels %d elapsed %.3f\n', ...
   k,min(labels(:,8)),size(labels,1),toc(timer)); end
end
[cost,index]=min(labels(:,8)); arrivals=zeros(n,1); departures=arrivals; v=zeros(n,3);
if ~isempty(incumbent)
 assert(cost<=incumbent.nominal_dv_km_s+1e-8,'zhangdp:incumbentCost', ...
  'Protected incumbent must remain available in the delayed graph.');
end
for k=n:-1:1
 row=layers{k}(index,:); arrivals(k)=row(1); departures(k)=row(10); v(k,:)=row(11:13); index=row(9);
end
path=struct('nominal_dv_km_s',cost,'arrival_times_s',arrivals,'departure_times_s',departures, ...
 'v_depart',v,'kind','fixed_order_delayed_finite_beam_two_body_proposal');
stats.elapsed_s=toc(timer);
end
