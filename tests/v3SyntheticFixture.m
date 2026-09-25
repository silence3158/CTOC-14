function [eph,c,s]=v3SyntheticFixture(mode)
%V3SYNTHETICFIXTURE 35 co-orbital synthetic targets; NEVER competition data.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
file=fullfile(sim,'runs','v3','preprocessing','targets_20260923_release','target_ephemeris.mat');
eph=ctocscreen.v3LoadTargetEphemeris(file); eph.synthetic=true;
eph.model.horizon_s=600; c=ctocscreen.v3Defaults(struct('scan_step_s',60,'max_step_s',30, ...
 'shooting_reltol',1e-12,'joint_iterations',3,'joint_seconds',10));
q=[eph.model.re+600 0.0001 -0.0002 1.1 .3 .4]; x=ctocscreen.initialState(q,eph.model.mu,eph.model.re);
eph.states0=repmat(x,35,1);
[~,~,sol]=ctocscreen.v3Arc(x,0,600,eph.model,c,false,true);
h=10; knots=(0:h:600).'; y=deval(sol,knots).'; a=ctocscreen.v3J2Acceleration(knots,y(:,1:3),eph.model);
r0=y(1:end-1,1:3); r1=y(2:end,1:3); v0=h*y(1:end-1,4:6); v1=h*y(2:end,4:6);
a0=h*h*a(1:end-1,:); a1=h*h*a(2:end,:); D=r1-r0-v0-a0/2; E=v1-v0-a0; F=a1-a0;
coef=cat(3,r0,v0,a0/2,10*D-4*E+F/2,-15*D+7*E-F,6*D-3*E+F/2);
eph.targets=repmat({struct('step_s',h,'coef',coef)},1,35);
s=struct('schema_version','free_maneuver_v3','dynamics_id','central_j2','initial_q',q, ...
 'maneuver_times_s',zeros(0,1),'delta_v_km_s',zeros(0,3),'duration_s',600, ...
 'witness_times_s',300*ones(35,1),'validation_level','proposal');
if nargin>0&&strcmp(mode,'separated_visits')
 s.witness_times_s=linspace(50,550,35).';
 opt=odeset('RelTol',3e-14,'AbsTol',1e-12,'MaxStep',30);
 for id=1:35
  t=s.witness_times_s(id); target=deval(sol,t);
  target(4:6)=target(4:6)+.05*[cos(id);sin(id);.5];
  back=ode89(@rhs,[t 0],target,opt); eph.states0(id,:)=deval(back,0).';
  forward=ode89(@rhs,[0 600],eph.states0(id,:).',opt); yy=deval(forward,knots).';
  aa=ctocscreen.v3J2Acceleration(knots,yy(:,1:3),eph.model);
  r0=yy(1:end-1,1:3); r1=yy(2:end,1:3); v0=h*yy(1:end-1,4:6); v1=h*yy(2:end,4:6);
  a0=h*h*aa(1:end-1,:); a1=h*h*aa(2:end,:); D=r1-r0-v0-a0/2; E=v1-v0-a0; F=a1-a0;
  eph.targets{id}.coef=cat(3,r0,v0,a0/2,10*D-4*E+F/2,-15*D+7*E-F,6*D-3*E+F/2);
 end
end
 function dz=rhs(t,z)
  dz=[z(4:6);ctocscreen.v3ReferenceForce(t,z(1:3).',eph.model).'];
 end
end
