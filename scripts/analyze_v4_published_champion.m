function result=analyze_v4_published_champion(label)
%ANALYZE_V4_PUBLISHED_CHAMPION Offline file diagnosis, never search input.
% Infer departure impulses by backward natural propagation of cached final
% states. Position closure tests the interpretation; fixed-control replay
% separately tests the inferred schedule without resetting any states.
if nargin<1, label='champion_20260928_01'; end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
folder=fullfile(sim,'runs/v4/diagnostics',label); assert(~isfolder(folder)); mkdir(folder);
source=fullfile(sim,'乙组前八名结果文件','乙组前八名结果文件');
files=dir(fullfile(source,'*.atk')); [~,ix]=sort({files.name}); files=files(ix);
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v4.defaults(); signature=ctocscreen.v4.signature(); clock=tic;
xp=javax.xml.xpath.XPathFactory.newInstance().newXPath();
result=struct('purpose','published_file_diagnosis_only_no_search_no_ATK', ...
 'signature',signature,'files',{{}},'source_directory',source,'official_engine_verified',false);
for fi=1:numel(files)
 file=fullfile(source,files(fi).name); doc=xmlread(file);
 seq=one(doc,'//Satellite[@Name="Inspector"]/Orbit/Propagator/Planning/Segment');
 ns=xp.evaluate('Segment',seq,javax.xml.xpath.XPathConstants.NODESET); N=ns.getLength();
 rec=struct('file',files(fi).name,'sha256',ctocscreen.v3FileHash(file), ...
  'segment_names',{{}},'segment_types',{{}},'duration_s',zeros(N,1), ...
  'cached_dv_km_s',nan(N,1),'target_ids',zeros(N,1),'method',nan(N,1));
 for k=1:N
  n=ns.item(k-1); rec.segment_names{k}=char(n.getAttribute('Name')); rec.segment_types{k}=char(n.getAttribute('ComponentType'));
  v=optional(n,'Duration'); if ~isempty(v), rec.duration_s(k)=str2double(v); end
  cp=xp.evaluate('CoordSystem/CoordAxes/CoordPoint',n,javax.xml.xpath.XPathConstants.NODE);
  if ~isempty(cp), id=sscanf(char(cp.getAttribute('Name')),'Target%d'); if ~isempty(id), rec.target_ids(k)=id; end, end
  dv=xp.evaluate('DV1',n,javax.xml.xpath.XPathConstants.NODE);
  if ~isempty(dv), rec.cached_dv_km_s(k)=str2double(char(dv.getAttribute('DV')))/1000; end
  v=optional(n,'MethodSwitch'); if ~isempty(v), rec.method(k)=str2double(v); end
 end
 [rec.unique_types,~,it]=unique(rec.segment_types); rec.type_counts=accumarray(it(:),1);
 result.files{fi}=rec;
 fprintf('FILE %d %s types=%s counts=%s cachedDVsum=%.9f\n',fi,rec.file,strjoin(rec.unique_types,','),mat2str(rec.type_counts.'),sum(rec.cached_dv_km_s,'omitnan'));
 if fi~=1, continue; end
 assert(strcmp(rec.segment_types{1},'CMCSInitialState')&&all(strcmp(rec.segment_types(2:end),'CMCSLambertTarget')));
 S=N-1; xb=zeros(6,S); xe=xb; times=[0;cumsum(rec.duration_s(2:end))]; offsets=zeros(S,1);
 prePostContinuity=zeros(S,2); cachedUTCResidual=zeros(S,2); targetDelta=zeros(35,1);
 epoch=datetime('2035-01-01 12:00:00','InputFormat','yyyy-MM-dd HH:mm:ss');
 x0=state(one(ns.item(0),'InitialState')); prior=x0;
 for k=1:S
  n=ns.item(k); assert(strcmp(optional(n,'Active'),'1')&&strcmp(optional(n,'IsOnedv'),'1'));
  xb(:,k)=state(one(n,'InitialState')); xe(:,k)=state(one(n,'FinalState'));
  prePostContinuity(k,:)=[norm(xb(1:3,k)-prior(1:3)),norm(xb(4:6,k)-prior(4:6))]; prior=xe(:,k);
  pos=one(n,'DesiredState/Pos'); offsets(k)=norm(arrayfun(@(z)str2double(char(pos.getAttribute(z))),{'X','Y','Z'}))/1000;
  for side=1:2
   paths={'InitialState/UTC','FinalState/UTC'};
   at=datetime(optional(n,paths{side}),'InputFormat','yyyy-MM-dd HH:mm:ss.SSSSSS');
   cachedUTCResidual(k,side)=seconds(at-epoch)-times(k+side-1);
  end
 end
 for id=1:35
  o=one(doc,sprintf('//Satellite[@Name="Target%d"]/Orbit',id));
  names={'PositionX','PositionY','PositionZ','VelocityX','VelocityY','VelocityZ'};
  vals=cellfun(@(z)str2double(optional(o,z))/1000,names);
  targetDelta(id)=max(abs(vals-eph.states0(id,:)));
 end
 result.champion=struct('x0',x0,'before',xb,'cached_end',xe,'times',times, ...
  'target_ids',rec.target_ids(2:end),'names',{rec.segment_names(2:end)}, ...
  'desired_offset_km',offsets,'continuity',prePostContinuity,'utc_residual_s',cachedUTCResidual, ...
  'target_state_max_difference',max(targetDelta),'initial_orbit',ctocscreen.v4.initialOrbit(x0,eph.model));
 ids=result.champion.target_ids; goalDistances=nan(S,1);
 for k=find(ids>0).'
  rt=ctocscreen.v3QueryTargets(eph,ids(k),times(k+1),'pairs');
  goalDistances(k)=norm(xe(1:3,k)-rt(:));
 end
 result.champion.cached_target_distances_km=goalDistances;
 modes={'nominal_pole','fixed_z','two_body'}; recon=cell(1,3);
 opt=odeset('RelTol',3e-14,'AbsTol',[1e-12*ones(3,1);1e-15*ones(3,1)],'MaxStep',120);
 for mi=1:3
  m=eph.model;
  if mi==2, m.pole=repmat([0 0 1],size(m.pole,1),1); end
  if mi==3, m.j2=0; end
  u=zeros(S,3); rerr=zeros(S,1); post=zeros(6,S);
  for k=1:S
   sol=ode89(@(t,x)[x(4:6);ctocscreen.v3ReferenceForce(t,x(1:3).',m).'],[times(k+1),times(k)],xe(:,k),opt);
   z=deval(sol,times(k)); post(:,k)=z;
   rerr(k)=norm(z(1:3)-xb(1:3,k)); u(k,:)=(z(4:6)-xb(4:6,k)).';
  end
  recon{mi}=struct('mode',modes{mi},'u',u,'post',post,'position_closure_km',rerr, ...
   'J',sum(vecnorm(u,2,2)),'max_closure',max(rerr),'median_closure',median(rerr));
  fprintf('REVERSE %s J=%.12f closure max=%.9g median=%.9g counts >1m/s=%d >10m/s=%d\n', ...
   modes{mi},recon{mi}.J,max(rerr),median(rerr),sum(vecnorm(u,2,2)>.001),sum(vecnorm(u,2,2)>.01));
 end
 result.champion.reconstructions=recon;
 save(fullfile(folder,'analysis.mat'),'result','-v7.3');
end
result.elapsed_s=toc(clock); result.source_unchanged=isequal(signature,ctocscreen.v4.signature()); assert(result.source_unchanged);
save(fullfile(folder,'analysis.mat'),'result','-v7.3');
fprintf('SAVED %s elapsed %.3f\n',folder,result.elapsed_s);

 function foundNode=one(context,path)
  foundNode=xp.evaluate(path,context,javax.xml.xpath.XPathConstants.NODE); assert(~isempty(foundNode),'Missing XML %s',path);
 end
 function s=optional(context,path)
  optionalNode=xp.evaluate(path,context,javax.xml.xpath.XPathConstants.NODE); s=''; if ~isempty(optionalNode), s=char(optionalNode.getTextContent()); end
 end
 function stateVector=state(stateNode)
  posNode=one(stateNode,'Pos'); velNode=one(stateNode,'Vel'); axes={'X','Y','Z'}; stateVector=zeros(6,1);
  for axisIndex=1:3, stateVector(axisIndex)=str2double(char(posNode.getAttribute(axes{axisIndex})))/1000; stateVector(axisIndex+3)=str2double(char(velNode.getAttribute(axes{axisIndex})))/1000; end
 end
end
