function [s,name]=v3Mutate(s,eph,c,stream,operator)
%V3MUTATE Wide changes affect the whole candidate, then require joint repair.
scale=1;
if nargin<5
 operator=randi(stream,8);
 if rand(stream)>c.wide_mutation_probability, scale=c.mutation_scale; end
end
T=s.duration_s; M=numel(s.maneuver_times_s); names={'insert','delete','move','initial_orbit', ...
 'reassign_witness','coast_duration','rebuild_window','merge'}; name=names{operator};
switch operator
 case 1
  if M<c.max_maneuvers
   inserted=rand(stream)*T;
   preceding=find(s.maneuver_times_s<inserted,1,'last');
   if ~isempty(preceding)
    s.delta_v_km_s(preceding,:)=s.delta_v_km_s(preceding,:)+randn(stream,1,3)*.2*scale;
   end
   s.maneuver_times_s(end+1,1)=inserted;
   s.delta_v_km_s(end+1,:)=randn(stream,1,3)*.2*scale;
  end
 case 2
  if M>0
   j=randi(stream,M);
   s.maneuver_times_s(j)=[]; s.delta_v_km_s(j,:)=[];
  end
 case 3
  if M>0, j=randi(stream,M); s.maneuver_times_s(j)=(1-scale)*s.maneuver_times_s(j)+scale*rand(stream)*T; end
 case 4
  root=ctocscreen.v3Roots(eph,ctocscreen.v3Defaults(struct('root_count',1)),stream);
  q=root{1}.schedule.initial_q; delta=q-s.initial_q;
  delta(5:6)=atan2(sin(delta(5:6)),cos(delta(5:6)));
  s.initial_q=s.initial_q+scale*delta;
 case 5
  known=find(isfinite(s.witness_times_s));
  if numel(known)>=2
   [~,order]=sort(s.witness_times_s(known)); ordered=known(order);
   a=randi(stream,numel(ordered)-1); pair=ordered(a:a+1);
   if scale==1, pair=known(randperm(stream,numel(known),2)); end
   s.witness_times_s(pair)=flipud(s.witness_times_s(pair));
   if s.witness_times_s(pair(1))==s.witness_times_s(pair(2))
    s.witness_times_s(pair)=sort(rand(stream,2,1)*T);
   end
  end
 case 6
  newT=min(eph.model.horizon_s,max(1,T*exp(randn(stream)*.5*scale)));
  unknown=isnan(s.witness_times_s);
  s.maneuver_times_s=min(newT,max(0,s.maneuver_times_s*(newT/T)));
  s.witness_times_s=min(newT,max(0,s.witness_times_s*(newT/T)));
  s.witness_times_s(unknown)=NaN;
  s.duration_s=newT;
 case 7
  if scale<1
   s.delta_v_km_s=s.delta_v_km_s+scale*.3*randn(stream,M,3);
   s.maneuver_times_s=(1-scale)*s.maneuver_times_s+scale*rand(stream,M,1)*T;
  else
  ab=sort(rand(stream,2,1)*T); keep=s.maneuver_times_s<ab(1)|s.maneuver_times_s>ab(2);
  s.maneuver_times_s=s.maneuver_times_s(keep); s.delta_v_km_s=s.delta_v_km_s(keep,:);
  count=min(randi(stream,4),max(0,c.max_maneuvers-sum(keep)));
  s.maneuver_times_s=[s.maneuver_times_s;ab(1)+diff(ab)*rand(stream,count,1)];
  s.delta_v_km_s=[s.delta_v_km_s;randn(stream,count,3)*.3];
  change=s.witness_times_s>=ab(1)&s.witness_times_s<=ab(2);
  s.witness_times_s(change)=ab(1)+diff(ab)*rand(stream,sum(change),1);
  end
 case 8
  if M>=2
   [~,ii]=sort(s.maneuver_times_s); j=randi(stream,M-1); a=ii(j); b=ii(j+1);
   s.maneuver_times_s(a)=(1-scale)*s.maneuver_times_s(a)+scale*mean(s.maneuver_times_s([a b]));
   s.delta_v_km_s(a,:)=s.delta_v_km_s(a,:)+s.delta_v_km_s(b,:);
   s.maneuver_times_s(b)=[]; s.delta_v_km_s(b,:)=[];
  end
end
if isfield(s,'visit_plan_times_s'), s.visit_plan_times_s=s.witness_times_s; end
s.validation_level='proposal'; s.mutation_scale=scale; s=ctocscreen.v3Normalize(s,eph.model);
end
