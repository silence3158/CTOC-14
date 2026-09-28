function out=analyze_v4_scene_controls(label)
%ANALYZE_V4_SCENE_CONTROLS Inspect published control knots and intervisit legs.
% Uses archived independent replay; no propagation or optimization.
if nargin<1, label='scene_controls_20260928_02'; end
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
p=load(fullfile(sim,'runs/v4/diagnostics/champion_gap_20260928_02/gap_analysis.mat'));
r=p.report; folder=fullfile(sim,'runs/v4/diagnostics',label);
assert(~isfolder(folder)); mkdir(folder);
xp=javax.xml.xpath.XPathFactory.newInstance().newXPath(); out=struct('ranks',{{}});
for rank=1:2
 d=r.trajectories{rank}; file=fullfile(sim,'乙组前八名结果文件','乙组前八名结果文件',r.ranked_files{rank}.file);
 assert(strcmp(ctocscreen.v3FileHash(file),r.ranked_files{rank}.sha256));
 doc=xmlread(file); seq=xp.evaluate('//Satellite[@Name="Inspector"]/Orbit/Propagator/Planning/Segment',doc,javax.xml.xpath.XPathConstants.NODE);
 ns=xp.evaluate('Segment',seq,javax.xml.xpath.XPathConstants.NODESET);
 rows=struct('name',{},'time_s',{},'duration_to_next_s',{},'dv_m_s',{},'kind',{},'end_label',{}); t=0;
 for j=1:ns.getLength()-1
  node=ns.item(j); kind=char(node.getAttribute('ComponentType')); name=char(node.getAttribute('Name'));
  if rank==1
   assert(strcmp(kind,'CMCSLambertTarget')); dt=number(node,'Duration');
   rows(end+1)=struct('name',name,'time_s',t,'duration_to_next_s',dt, ...
    'dv_m_s',1000*norm(d.q.u(j,:)),'kind',kind,'end_label',name); %#ok<AGROW>
   t=t+dt;
  elseif strcmp(kind,'CMCSManeuver')
   burn=xp.evaluate('Burn',node,javax.xml.xpath.XPathConstants.NODE); axes={'X','Y','Z'};
   dv=cellfun(@(a)str2double(char(burn.getAttribute(a))),axes);
   rows(end+1)=struct('name',name,'time_s',t,'duration_to_next_s',NaN, ...
    'dv_m_s',norm(dv),'kind',kind,'end_label',''); %#ok<AGROW>
  elseif strcmp(kind,'CMCSPropagate')
   trigger=xp.evaluate('StopCondition/StateTrigger',node,javax.xml.xpath.XPathConstants.NODE);
   dt=str2double(char(trigger.getAttribute('TripValue'))); t=t+dt;
   rows(end).duration_to_next_s=dt; rows(end).end_label=name;
  else
   error('Unsupported scene segment %s',kind);
  end
 end
 raw=struct2table(rows); assert(abs(t-d.q.T)<1e-5);
 wt=d.visits.day*86400; ids=d.visits.target;
 raw.nearest_visit_s=min(abs(raw.time_s-wt.'),[],2);
 raw.at_visit_60s=raw.nearest_visit_s<=60; raw.major=raw.dv_m_s>1;
 raw.hours=raw.time_s/3600;
 counts=zeros(35,1); inside=counts; visitBurn=counts; dvlegs=counts;
 from=[0;ids(1:end-1)]; starts=[0;wt(1:end-1)];
 lists=cell(35,1); interiorLists=lists;
 for k=1:35
  % A pulse within +/-60s of a visit belongs to the departure AFTER that visit.
  lo=starts(k)-60; if k==1, lo=-Inf; end
  hi=wt(k)-60;
  ix=find(raw.time_s>=lo&raw.time_s<hi); im=ix(raw.major(ix));
  ii=im(raw.time_s(im)>starts(k)+60);
  if k==1, ii=im; end
  counts(k)=numel(im); inside(k)=numel(ii); visitBurn(k)=counts(k)-inside(k);
  dvlegs(k)=sum(raw.dv_m_s(ix)); lists{k}=raw.name(im); interiorLists{k}=raw.name(ii);
 end
 legs=table(from,ids,starts/3600,wt/3600,(wt-starts)/3600,counts,inside,visitBurn,dvlegs,lists,interiorLists, ...
  'VariableNames',{'from_target','to_target','start_h','end_h','duration_h','major_count','interior_major_count','departure_visit_major_count','dv_m_s','major_names','interior_names'});
 norms=sort(raw.dv_m_s,'descend'); caps=[1 5 10 20 30 40];
 cumulative=arrayfun(@(k)sum(norms(1:min(k,numel(norms))))/sum(norms),caps);
 stats=struct('raw_count',height(raw),'major_count',sum(raw.major),'total_dv_m_s',sum(raw.dv_m_s), ...
  'small_count',sum(~raw.major),'small_total_m_s',sum(raw.dv_m_s(~raw.major)), ...
  'major_interval_quantiles_h',prctile(diff(raw.hours(raw.major)),[0 25 50 75 100]), ...
  'raw_interval_quantiles_h',prctile(diff(raw.hours),[0 25 50 75 100]), ...
  'raw_interval_cv',std(diff(raw.hours))/mean(diff(raw.hours)), ...
  'major_near_visit_count',sum(raw.major&raw.at_visit_60s), ...
  'major_far_visit_count',sum(raw.major&~raw.at_visit_60s), ...
  'intervisit_major_histogram_0_to_max',arrayfun(@(k)sum(counts(2:end)==k),0:max(counts(2:end))), ...
  'intervisit_interior_major_histogram_0_to_max',arrayfun(@(k)sum(inside(2:end)==k),0:max(inside(2:end))), ...
  'top_k',caps,'top_k_cost_fraction',cumulative);
 x=d.q.x0; h=cross(x(1:3),x(4:6));
 stats.initial_inc_deg=acosd(h(3)/norm(h)); stats.initial_raan_deg=mod(atan2d(h(1),-h(2)),360);
 stats.distance_quantiles_km=prctile(d.verification.distance_km,[0 25 50 75 100]);
 bm=d.burns(d.burns.dv_km_s>.001,:);
 stats.major_normal_fraction_quantiles=prctile(abs(bm.normal)./bm.dv_km_s,[0 25 50 75 100]);
 stats.major_tangential_fraction_quantiles=prctile(abs(bm.tangential)./bm.dv_km_s,[0 25 50 75 100]);
 stats.major_dv_after_initial_km_s=sum(bm.dv_km_s(2:end));
 out.ranks{rank}=struct('file',file,'sha256',r.ranked_files{rank}.sha256,'raw',raw,'legs',legs,'stats',stats);
 writetable(raw,fullfile(folder,sprintf('rank%d_scene_controls.csv',rank)));
 writetable(removevars(legs,{'major_names','interior_names'}),fullfile(folder,sprintf('rank%d_intervisit_counts.csv',rank)));
 fprintf('RANK %d\n',rank); disp(stats); disp(legs(:,1:9));
end
save(fullfile(folder,'scene_controls.mat'),'out');
fid=fopen(fullfile(folder,'scene_controls.json'),'w','n','UTF-8'); cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(out));

 function x=number(context,path)
  selected=xp.evaluate(path,context,javax.xml.xpath.XPathConstants.NODE); assert(~isempty(selected));
  x=str2double(char(selected.getTextContent()));
 end
end
