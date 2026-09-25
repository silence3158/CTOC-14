function roots=v3Roots(eph,c,stream)
%V3ROOTS Actual six-dimensional initial-orbit diversity, not RNG labels alone.
roots=cell(1,c.root_count);
for k=1:c.root_count
 e=.000999*sqrt(rand(stream)); w=2*pi*rand(stream);
 q=[eph.model.re+590+20*rand(stream),e*cos(w),e*sin(w), ...
  acos(2*rand(stream)-1),2*pi*rand(stream),2*pi*rand(stream)];
 target=0; duration=NaN; method='random'; wait=0;
 if mod(k,2)==0, wait=min(.25*eph.model.horizon_s,c.root_wait_max_s)*rand(stream); end
 if rand(stream)<c.guided_root_fraction
  target=randi(stream,35); r0=ctocscreen.v3QueryTargets(eph,target,0);
  duration=min(.8*eph.model.horizon_s,pi*sqrt(((q(1)+norm(r0))/2)^3/eph.model.mu));
  duration=max(1,min(eph.model.horizon_s-wait,duration*(.9+.2*rand(stream))));
  [r,v]=ctocscreen.v3QueryTargets(eph,target,wait+duration); h=cross(r,v); h=h/norm(h);
  inc=acos(max(-1,min(1,h(3)))); Omega=atan2(h(1),-h(2));
  p=[cos(Omega),sin(Omega),0]; b=cross(h,p);
  launch=-r/norm(r); phase=(2*rand(stream)-1)*.15;
  q(4:6)=[inc,Omega,atan2(dot(launch,b),dot(launch,p))+phase];
  method='target_geometry';
 end
 if c.cost_guidance_enabled&&k<=ceil(c.root_count/2)
  % Low-inclination starts are a proposal family, not a fixed initial plane.
  radii=vecnorm(eph.states0(:,1:3),2,2);
  angular=cross(eph.states0(:,1:3),eph.states0(:,4:6),2);
  inclinations=acosd(angular(:,3)./vecnorm(angular,2,2));
  low=find(inclinations<5); if isempty(low), low=(1:35).'; end
  [~,sorted]=sort(radii(low)); near=low(sorted);
  target=near(1+mod(floor((k-1)/2),min(2,numel(near))));
  duration=pi*sqrt(((q(1)+radii(target))/2)^3/eph.model.mu);
  duration=min([.8*eph.model.horizon_s,eph.model.horizon_s-wait,duration]);
  r=ctocscreen.v3QueryTargets(eph,target,wait+duration); direction=r/norm(r);
  minimumTilt=asin(abs(direction(3)));
  tilt=deg2rad(c.root_plane_tilts_deg(1+mod(k-1,numel(c.root_plane_tilts_deg))))+deg2rad(2*rand(stream));
  tilt=max(minimumTilt+1e-8,min(pi-minimumTilt-1e-8,tilt));
  pole=[0 0 1]-direction(3)*direction; pole=pole/norm(pole);
  turn=acos(max(-1,min(1,cos(tilt)/pole(3))));
  if rand(stream)<.5, turn=-turn; end
  h=cos(turn)*pole+sin(turn)*cross(direction,pole);
  q(4)=acos(max(-1,min(1,h(3)))); q(5)=atan2(h(1),-h(2));
  p=[cos(q(5)),sin(q(5)),0]; b=cross(h,p);
  launch=-direction;
  q(6)=atan2(dot(launch,b),dot(launch,p))+.03*(2*rand(stream)-1);
  method='plane_diverse_geometry';
 end
 % Secular rates only transport a guess. Waiting is always replayed in J2.
 n=sqrt(eph.model.mu/q(1)^3); semilatus=q(1)*(1-e^2);
 scale=eph.model.j2*(eph.model.re/semilatus)^2; ci=cos(q(4));
 nodeRate=-1.5*scale*n*ci;
 latitudeRate=n+.75*scale*n*((5*ci^2-1)+sqrt(1-e^2)*(3*ci^2-1));
 q(5)=q(5)-nodeRate*wait; q(6)=q(6)-latitudeRate*wait;
 key=sprintf('%.17g/',q);
 s=struct('schema_version','free_maneuver_v3','dynamics_id','central_j2', ...
  'initial_q',q,'maneuver_times_s',zeros(0,1),'delta_v_km_s',zeros(0,3), ...
  'duration_s',1,'witness_times_s',nan(35,1),'validation_level','proposal','root_key',key);
 roots{k}=struct('schedule',s,'t',0,'state',ctocscreen.initialState(q,eph.model.mu,eph.model.re), ...
  'visited',false(35,1),'J',0,'keys',{{}},'last_gain',0,'estimate',0, ...
  'seed_target',target,'seed_duration_s',duration,'seed_wait_s',wait,'root_method',method);
end
end
