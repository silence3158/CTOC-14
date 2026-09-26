function node=root(index,eph,c,stream)
%ROOT Fresh physical initial conditions from raw target geometry or free draws.
m=eph.model; a=m.re+590+20*rand(stream); e=.00095*sqrt(rand(stream)); w=2*pi*rand(stream);
q0=[a,e*cos(w),e*sin(w),acos(2*rand(stream)-1),2*pi*rand(stream),2*pi*rand(stream)];
radii=vecnorm(eph.states0(:,1:3),2,2); [~,rank]=sort(radii);
target=rank(1+mod(index-1,numel(rank))); duration=pi*sqrt(((a+radii(target))/2)^3/m.mu);
kind='free';
if mod(index,4)~=0
 if mod(index,3)==2, target=1+mod(7*index-1,35); end
 duration=pi*sqrt(((a+radii(target))/2)^3/m.mu);
 [r,v]=ctocscreen.v3QueryTargets(eph,target,duration);
 h=cross(r,v); h=h/norm(h); inc=acos(max(-1,min(1,h(3)))); om=atan2(h(1),-h(2));
 p=[cos(om),sin(om),0]; b=cross(h,p); direction=-r/norm(r);
 q0(4:6)=[inc,om,atan2(dot(direction,b),dot(direction,p))+.03*(2*rand(stream)-1)];
 kind='raw_target_geometry';
end
q=struct('x0',ctocscreen.initialState(q0,m.mu,m.re).','tau',zeros(0,1), ...
 'u',zeros(0,3),'T',0,'witness',nan(35,1));
[actual,trace]=ctocscreen.v4.replay(q,eph,c);
node=struct('q',q,'actual',actual,'trace',trace,'root_id',index, ...
 'root_kind',kind,'seed_target',target,'seed_duration',duration,'attempts',0, ...
 'zero_gain',0,'generation',0,'origin','root','heuristic_H',NaN);
end
