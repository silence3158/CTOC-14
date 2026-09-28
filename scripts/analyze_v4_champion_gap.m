function report=analyze_v4_champion_gap(label)
%ANALYZE_V4_CHAMPION_GAP Fixed-control diagnostics of published trajectories.
% No optimization, re-aiming, search initialization or original file changes.
if nargin<1, label='champion_gap_20260928_01'; end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
folder=fullfile(sim,'runs/v4/diagnostics',label); assert(~isfolder(folder)); mkdir(folder);
raw=load(fullfile(sim,'runs/v4/diagnostics/champion_20260928_02/analysis.mat')); raw=raw.result;
source=raw.source_directory; xp=javax.xml.xpath.XPathFactory.newInstance().newXPath();
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v4.defaults(); signature=ctocscreen.v4.signature(); clock=tic;
report=struct('purpose','postcompetition_fixed_control_diagnosis_not_search', ...
 'signature',signature,'champion_extraction',raw.champion,'ranked_files',{raw.files},'trajectories',{{}});
for fi=1:numel(raw.files)
 rec=raw.files{fi}; file=fullfile(source,rec.file); assert(strcmp(rec.sha256,ctocscreen.v3FileHash(file)));
 if ~any(strcmp(rec.segment_types,'CMCSManeuver')), continue; end
 doc=xmlread(file); seq=one(doc,'//Satellite[@Name="Inspector"]/Orbit/Propagator/Planning/Segment');
 ns=xp.evaluate('Segment',seq,javax.xml.xpath.XPathConstants.NODESET);
 burnNorm=[]; burnTime=[]; t=0; impulses=zeros(0,3);
 for k=0:ns.getLength()-1
  node=ns.item(k); kind=char(node.getAttribute('ComponentType'));
  assert(strcmp(textAt(node,'Active'),'1'));
  if strcmp(kind,'CMCSManeuver')
   settings=one(node,'OrbBurn'); assert(strcmp(char(settings.getAttribute('BurnType')),'0'));
   assert(strcmp(char(settings.getAttribute('BurnCoordAxes')),'ICRF')&&strcmp(char(settings.getAttribute('BurnCoordType')),'0'));
   dv=xyz(one(node,'Burn'))/1000; burnNorm(end+1,1)=norm(dv); burnTime(end+1,1)=t; impulses(end+1,:)=dv.'; %#ok<AGROW>
  elseif strcmp(kind,'CMCSPropagate')
   trigger=one(node,'StopCondition/StateTrigger'); t=t+str2double(char(trigger.getAttribute('TripValue')));
  elseif strcmp(kind,'CMCSInitialState')
   assert(k==0); st=one(node,'InitialState'); xstart=[xyz(one(st,'Pos'));xyz(one(st,'Vel'))]/1000;
  else
   error('Unrecognized segment %s',kind);
  end
 end
 rec.explicit_burn_count=numel(burnNorm); rec.explicit_J=sum(burnNorm); rec.duration_from_controls=t;
 rec.sub_1mps_count=sum(burnNorm<.001); rec.major_10mps_count=sum(burnNorm>.01);
 rec.maximum_burn=max(burnNorm); report.ranked_files{fi}=rec;
 if fi==2
  runnerq=struct('x0',xstart,'tau',burnTime,'u',impulses,'T',t,'witness',nan(35,1));
 end
end
z=raw.champion; q=struct('x0',z.x0,'tau',z.times(1:end-1), ...
 'u',z.reconstructions{1}.u,'T',z.times(end),'witness',nan(35,1));
for k=find(z.target_ids>0).', q.witness(z.target_ids(k))=z.times(k+1); end
warm=load(fullfile(sim,'runs/v4/development/incumbent_shared_warm_412_300_seed888_01/report.mat'));
ours=load(fullfile(sim,'runs/v4/development/incumbent_breadth_1800_seed888_01/report.mat'));
inputs={q,runnerq,warm.report.final.q,ours.report.final.q};
labels={'rank1_cached_state_inferred','rank2_explicit_controls','teammate_reconstructed','our_cold_1800'};
for a=1:4
 qt=inputs{a}; previous=fullfile(sim,'runs/v4/diagnostics/champion_gap_20260928_01/gap_analysis.mat');
 if a==1&&isfile(previous)
  prior=load(previous,'report'); d=prior.report.trajectories{1};
  assert(isequal(qt.x0,d.q.x0)&&isequal(qt.tau,d.q.tau)&&isequal(qt.u,d.q.u)&&qt.T==d.q.T);
  v=d.verification; d.reused_identical_diagnostic=previous;
 else
  rt=tic; [v,tr]=ctocscreen.v4.replay(qt,eph,c,true);
  d=struct('label',labels{a},'q',tr.q,'input_q',qt,'verification',v,'replay_seconds',toc(rt));
  assert(~strcmp(v.status,'propagation_failure'),'Replay error: %s',v.failure_reason);
  d=describe(d,tr,eph.model);
 end
 report.trajectories{a}=d;
 fprintf('REPLAY %s visits=%d J=%.12f height=%.6f maxdistance=%.6f time=%.3f burns=%d\n', ...
  labels{a},v.visit_count,v.total_dv_km_s,v.min_altitude_lower_km,max(v.distance_km),d.replay_seconds,numel(qt.tau));
 fprintf('STRUCTURE thresholds=%s counts=%s offvisitMajor=%d/%d zeroVisitMajorArcs=%d multiVisitMajorArcs=%d\n', ...
  mat2str(d.thresholds),mat2str(d.counts),sum(d.burns.nearest_visit_s(d.burns.dv_km_s>.001)>60),sum(d.burns.dv_km_s>.001), ...
  sum(d.major_arc_visit_counts==0),sum(d.major_arc_visit_counts>=2));
 disp(d.percentiles);
 writetable(d.burns,fullfile(folder,[labels{a} '_burns.csv']));
 writetable(d.visits,fullfile(folder,[labels{a} '_visits.csv']));
 save(fullfile(folder,'gap_analysis.mat'),'report','-v7.3');
end
report.elapsed_s=toc(clock); report.source_unchanged=isequal(signature,ctocscreen.v4.signature()); assert(report.source_unchanged);
save(fullfile(folder,'gap_analysis.mat'),'report','-v7.3');
small=report; small.champion_extraction=rmfield(small.champion_extraction,{'reconstructions','before','cached_end'});
for a=1:numel(small.trajectories), small.trajectories{a}=rmfield(small.trajectories{a},'q'); end
fid=fopen(fullfile(folder,'summary.json'),'w','n','UTF-8'); cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(small));
fprintf('SAVED %s elapsed %.3f\n',folder,report.elapsed_s);

 function found=one(context,path)
  found=xp.evaluate(path,context,javax.xml.xpath.XPathConstants.NODE); assert(~isempty(found),'Missing XML %s',path);
 end
 function out=textAt(context,path)
  selected=one(context,path); out=char(selected.getTextContent());
 end
 function out=xyz(context)
  out=zeros(3,1); axes={'X','Y','Z'};
  for j=1:3, out(j)=str2double(char(context.getAttribute(axes{j}))); end
 end
end

function d=describe(d,tr,m)
q=d.q; v=d.verification; N=numel(q.tau); dv=vecnorm(q.u,2,2);
alt=zeros(N,1); speed=alt; afterA=alt; afterE=alt; nearVisit=alt; radial=alt; normal=alt; tangential=alt;
for k=1:N
 ix=find(tr.times==q.tau(k),1); x=tr.pre(:,ix); y=tr.post(:,ix);
 R=x(1:3)/norm(x(1:3)); H=cross(x(1:3),x(4:6)); H=H/norm(H); T=cross(H,R);
 alt(k)=norm(x(1:3))-m.re; speed(k)=norm(x(4:6));
 afterA(k)=1/(2/norm(y(1:3))-dot(y(4:6),y(4:6))/m.mu);
 afterE(k)=norm(cross(y(4:6),cross(y(1:3),y(4:6)))/m.mu-y(1:3)/norm(y(1:3)));
 nearVisit(k)=min(abs(q.tau(k)-v.witness_times_s(v.distance_km<=1)));
 radial(k)=dot(q.u(k,:),R); normal(k)=dot(q.u(k,:),H); tangential(k)=dot(q.u(k,:),T);
end
d.thresholds=[1e-6 1e-4 .001 .01]; d.counts=arrayfun(@(z)sum(dv>z),d.thresholds);
d.burns=table((1:N).',q.tau/86400,dv,cumsum(dv),alt,speed,afterA,afterE,nearVisit,radial,tangential,normal, ...
 'VariableNames',{'index','day','dv_km_s','cumulative_dv','altitude_km','speed_km_s','after_a_km','after_e','nearest_visit_s','radial','tangential','normal'});
[wt,ord]=sort(v.witness_times_s); passed=v.distance_km(ord)<=1;
vcost=arrayfun(@(z)sum(dv(q.tau<=z)),wt);
d.visits=table((1:35).',ord,wt/86400,v.distance_km(ord),passed,vcost, ...
 'VariableNames',{'order','target','day','distance_km','passed','cumulative_dv'});
major=dv>.001; edges=unique([0;q.tau(major);q.T]); d.major_arc_bounds=edges;
validW=v.witness_times_s(v.distance_km<=1); counts=zeros(numel(edges)-1,1); sets=cell(size(counts));
for k=1:numel(counts)
 ids=find(v.distance_km<=1&v.witness_times_s>edges(k)&v.witness_times_s<=edges(k+1));
 sets{k}=ids.'; counts(k)=numel(ids);
end
d.major_arc_visit_counts=counts; d.major_arc_targets=sets;
d.first_burn_dv=dv(find(major,1)); d.major_dv_sum=sum(dv(major)); d.minor_dv_sum=sum(dv(~major));
d.percentiles=struct('major_burn_dv',prctile(dv(major),[0 25 50 75 100]), ...
 'major_burn_speed',prctile(speed(major),[0 25 50 75 100]), ...
 'major_burn_altitude',prctile(alt(major),[0 25 50 75 100]), ...
 'major_delay_to_nearest_visit_s',prctile(nearVisit(major),[0 25 50 75 100]));
d.cost_at_coverage_10_20_30_35=nan(1,4); levels=[10 20 30 35];
for k=1:4, if numel(validW)>=levels(k), ws=sort(validW); d.cost_at_coverage_10_20_30_35(k)=sum(dv(q.tau<=ws(levels(k)))); end, end
d.final_days=q.T/86400;
end
