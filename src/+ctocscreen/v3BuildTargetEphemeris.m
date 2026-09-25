function cache = v3BuildTargetEphemeris(folder,options)
%V3BUILDTARGETEPHEMERIS Build and validate all 35 targets over ten days.
% Quintic Hermite r/v/a endpoints; velocity = derivative of same polynomial.
% Independent ode89 reference uses a separately expressed force equation.
if nargin<2, options=struct(); end
defaults=struct('initial_step_s',120,'minimum_step_s',15, ...
 'position_tolerance_km',1e-4,'velocity_tolerance_km_s',1e-7, ...
 'build_reltol',3e-14,'reference_reltol',3e-14);
keys=fieldnames(options);
assert(all(ismember(keys,fieldnames(defaults))),'ctocscreen:v3:option','Unknown option.');
for j=1:numel(keys), defaults.(keys{j})=options.(keys{j}); end
options=defaults;
names=fieldnames(defaults);
for j=1:numel(names)
 validateattributes(options.(names{j}),{'double'},{'scalar','positive','finite','real'});
end
validateattributes(options.initial_step_s,{'double'},{'scalar','positive','finite'});
validateattributes(options.minimum_step_s,{'double'},{'scalar','positive','finite','<=',options.initial_step_s});
assert(mod(864000,options.initial_step_s)==0,'ctocscreen:v3:step','Initial step must divide horizon.');
assert(~isfolder(folder),'ctocscreen:v3:existingFolder','Choose a new run folder.');
mkdir(folder); started=tic;
[model,audit,p]=ctocscreen.v3TargetModel();
cache=struct('schema_version','target_ephemeris_v3_1','model',model,'audit',audit, ...
 'options',options,'target_ids',p.target_ids,'states0',p.states0, ...
 'units',struct('position','km','velocity','km/s','acceleration','km/s^2','time','seconds since epoch'), ...
 'signature',ctocscreen.v3PreprocessSignature(),'targets',{cell(1,35)}, ...
 'created_utc',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss')));
report=zeros(35,8);
for id=1:35
 timer=tic; x0=p.states0(id,:).';
 buildOpt=odeset('RelTol',options.build_reltol,'AbsTol',[1e-12*ones(3,1);1e-15*ones(3,1)],'MaxStep',60);
 refOpt=odeset('RelTol',options.reference_reltol,'AbsTol',[1e-12*ones(3,1);1e-15*ones(3,1)],'MaxStep',60);
 sol=ode113(@(t,x)buildRhs(t,x,model),[0 model.horizon_s],x0,buildOpt);
 ref=ode89(@(t,x)referenceRhs(t,x,model),[0 model.horizon_s],x0,refOpt);
 h=options.initial_step_s; attempts=[];
 while true
  knots=(0:h:model.horizon_s).'; y=deval(sol,knots).';
  acc=ctocscreen.v3J2Acceleration(knots,y(:,1:3),model);
  cache.targets{id}=struct('step_s',h,'coef',hermite(y,acc,h));
  % Three deterministic interior points in EVERY segment, not just random samples.
  tq=reshape(knots(1:end-1)+h*[0.211324865405187 0.5 0.788675134594813],[],1);
  tq=sort([tq;knots]);
  [r,v]=ctocscreen.v3QueryTargets(cache,id,tq);
  yr=deval(ref,tq).'; yb=deval(sol,tq).';
  er=vecnorm(r-yr(:,1:3),2,2); ev=vecnorm(v-yr(:,4:6),2,2);
  interpR=max(vecnorm(r-yb(:,1:3),2,2));
  solverR=max(vecnorm(yb(:,1:3)-yr(:,1:3),2,2));
  [maxR,ix]=max(er); maxV=max(ev);
  attempts(end+1,:)=[h,maxR,maxV,interpR,solverR]; %#ok<AGROW>
  if maxR<=options.position_tolerance_km && maxV<=options.velocity_tolerance_km_s
   break
  end
  assert(h/2>=options.minimum_step_s,'ctocscreen:v3:accuracy', ...
   'Target %d failed accuracy (%.6g km, %.6g km/s). No valid cache published.',id,maxR,maxV);
  h=h/2;
 end
 cache.targets{id}.validation=struct('max_position_error_km',maxR, ...
  'max_velocity_error_km_s',maxV,'max_interpolation_position_error_km',interpR, ...
  'max_integrator_position_difference_km',solverR,'worst_time_s',tq(ix), ...
  'sample_count',numel(tq),'attempts',attempts);
 report(id,:)=[id,h,maxR*1000,maxV*1000,interpR*1000,solverR*1000,numel(tq),toc(timer)];
 fprintf('Target %02d: step=%g s max error=%.6g m velocity=%.6g m/s (%.1fs)\n', ...
  id,h,maxR*1000,maxV*1000,report(id,8));
end
cache.validation=struct('passed',true,'level','nominal_j2_numerical_verified', ...
 'max_position_error_m',max(report(:,3)),'max_velocity_error_m_s',max(report(:,4)), ...
 'method','ode113 + quintic versus independent ode89; knots and 3 interior points per segment', ...
 'limitation','Sampled numerical comparison, not a rigorous uniform bound or ATK alignment', ...
 'build_seconds',toc(started));
cache.validation.table=array2table(report,'VariableNames',{'target_id','step_s','max_position_error_m', ...
 'max_velocity_error_m_s','interpolation_error_m','integrator_difference_m','samples','seconds'});
% Warm JIT and benchmark five representative public query modes.
ctocscreen.v3QueryTargets(cache,1,12345);
cache.benchmark=bench(cache);
assert(isequaln(cache.signature,ctocscreen.v3PreprocessSignature()), ...
 'ctocscreen:v3:sourceChanged','Sources changed during build; rebuild required.');
temporary=fullfile(folder,'target_ephemeris.pending.mat');
published=fullfile(folder,'target_ephemeris.mat');
save(temporary,'cache','-v7.3');
movefile(temporary,published);
checksum=ctocscreen.v3FileHash(published);
hashFile=fopen([published '.sha256'],'w');
assert(hashFile>=0,'ctocscreen:v3:write','Cannot write cache checksum.');
fprintf(hashFile,'%s\n',checksum); fclose(hashFile);
writetable(cache.validation.table,fullfile(folder,'validation.csv'));
fid=fopen(fullfile(folder,'report.md'),'w','n','UTF-8'); c=onCleanup(@()fclose(fid));
fprintf(fid,'# V3 target ephemeris preprocessing\n\nCreated UTC: %s\n\n',cache.created_utc);
fprintf(fid,'35 targets, [0,864000] s, km / s / km/s.\n\n');
fprintf(fid,'Frame: `%s`. Official ATK alignment: **pending**.\n\n',model.frame);
fprintf(fid,'Orientation: %s. Source EOP ends %s.\n\n',model.orientation,audit.eop_last_date);
fprintf(fid,'mu=%.12g km^3/s^2; RE=%.12g km; normalized C20=%.15g; J2=%.15g.\n\n',model.mu,model.re,model.c20_normalized,model.j2);
fprintf(fid,'CSV/source max component errors: position %.9g km, velocity %.9g km/s.\n\n', ...
 audit.max_position_component_error_km,audit.max_velocity_component_error_km_s);
fprintf(fid,'Max sampled state difference: %.9g m, %.9g m/s. Build %.3f s.\n\n', ...
 cache.validation.max_position_error_m,cache.validation.max_velocity_error_m_s,cache.validation.build_seconds);
fprintf(fid,'Reference: ode89, separately expressed central+J2 force, max step 60 s, RelTol %.3g; build ode113 RelTol %.3g.\n\n',options.reference_reltol,options.build_reltol);
fprintf(fid,'Pole table midpoint discrepancy: %.9g rad. Nominal pole model shared by both integrators.\n\n',model.pole_midpoint_error_rad);
fprintf(fid,'Not a rigorous uniform interpolation bound; no ATK process/library was executed.\n\n');
fprintf(fid,'Median query seconds (warmed, 100 repetitions per timing, 3 timings):\n\n');
fields=fieldnames(cache.benchmark);
for k=1:numel(fields), fprintf(fid,'- %s: %.9g\n',fields{k},cache.benchmark.(fields{k})); end
end

function dx=buildRhs(t,x,m)
dx=[x(4:6);ctocscreen.v3J2Acceleration(t,x(1:3).',m).'];
end

function dx=referenceRhs(t,x,m)
% Independently arranged radial/axial formula; separate integrator.
h=m.pole_times(2); i=min(floor(t/h)+1,size(m.pole,1)-1); u=(t-(i-1)*h)/h;
k=((1-u)*m.pole(i,:)+u*m.pole(i+1,:)).'; k=k/norm(k);
r=x(1:3); R=norm(r); axial=dot(r,k)*k; transverse=r-axial;
q=(dot(r,k)/R)^2; f=1.5*m.j2*(m.re/R)^2;
a=-m.mu/R^3*((1+f*(1-5*q))*transverse+(1+f*(3-5*q))*axial);
dx=[x(4:6);a];
end

function c=hermite(y,a,h)
r0=y(1:end-1,1:3); r1=y(2:end,1:3);
v0=h*y(1:end-1,4:6); v1=h*y(2:end,4:6);
a0=h^2*a(1:end-1,:); a1=h^2*a(2:end,:);
D=r1-r0-v0-a0/2; E=v1-v0-a0; F=a1-a0;
c=cat(3,r0,v0,a0/2,10*D-4*E+F/2,-15*D+7*E-F,6*D-3*E+F/2);
end

function b=bench(c)
f={@()ctocscreen.v3QueryTargets(c,1,12345), ...
 @()ctocscreen.v3QueryTargets(c,1,linspace(0,864000,1000)), ...
 @()ctocscreen.v3QueryTargets(c,[],12345), ...
 @()ctocscreen.v3QueryTargets(c,1:35,linspace(0,864000,35)), ...
 @()accQuery(c)};
names={'one_target','one_target_1000_times','all_targets_one_time','pairs_35','all_targets_acceleration'};
for k=1:numel(f)
 times=zeros(1,3);
 for j=1:3
  t=tic; for z=1:100, f{k}(); end; times(j)=toc(t)/100;
 end
 b.(names{k})=median(times);
end
end
function accQuery(c)
[~,~,~]=ctocscreen.v3QueryTargets(c,[],12345);
end
