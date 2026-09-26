function report=check_v4_crossing_t2(width,label,options)
%CHECK_V4_CROSSING_T2 Can the crossing trunk alone build a 35-target chain?
% Construction diagnostic, NOT a mission result: greedy layered beam of the
% given width over expand() only (no B, no shared arcs). Cold start, seed 888.
% Final chain is independently replayed so the report states actual visits.
if nargin<1, width=1; end
if nargin<2||isempty(label), label=char(datetime('now','Format','yyyyMMdd_HHmmss')); end
if nargin<3, options=struct(); end
options.beam_width=max(2,width);
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
eph=ctocscreen.v3LoadTargetEphemeris(fullfile(sim,'runs/v3/preprocessing/targets_20260923_release/target_ephemeris.mat'));
c=ctocscreen.v4.defaults(options); stream=RandStream('mt19937ar','Seed',888);
beam={}; for k=1:4, beam{end+1}=ctocscreen.v4.root(k,eph,c,stream); end %#ok<AGROW>
clock=tic; memory=[]; history=zeros(0,4); step=0;
while toc(clock)<900&&step<120
 step=step+1; next={};
 for k=1:min(width,numel(beam))
  [children,~,memory]=ctocscreen.v4.expand(beam{k},eph,c,stream,memory,c.action_seconds);
  next=[next,children]; %#ok<AGROW>
 end
 if isempty(next), break; end
 [pool,~]=ctocscreen.v4.selectBeam(next,c,eph);
 beam=pool(1:min(width,numel(pool)));
 b=beam{1}; history(end+1,:)=[toc(clock),b.actual.visit_count,b.actual.total_dv_km_s,b.q.T]; %#ok<AGROW>
 fprintf('step %3d t=%6.1f visits=%2d dv=%8.4f T=%.2f d H=%.2f\n',step,history(end,1),history(end,2),history(end,3),b.q.T/86400,b.heuristic_H);
 if b.actual.visit_count>=35||b.q.T>=eph.model.horizon_s-1, break; end
end
counts=cellfun(@(n)n.actual.visit_count,beam); [~,j]=max(counts); best=beam{j};
[verified,~]=ctocscreen.v4.replay(best.q,eph,c,true);
fprintf('T2 width=%d independent %d/35 dv=%.9f passed=%d height=%.3f km elapsed=%.1f s\n', ...
 width,verified.visit_count,verified.total_dv_km_s,verified.passed,verified.min_altitude_lower_km,toc(clock));
inc=verified.inclination_changes_deg; pl=verified.plane_angles_deg;
fprintf('plane rotation total=%.1f deg, max=%.1f, burns>5deg=%d of %d\n',sum(pl),max([pl;0]),sum(pl>5),numel(pl));
report=struct('diagnostic_only',true,'width',width,'options',options,'history',history,'verification',verified, ...
 'q',best.q,'signature',ctocscreen.v4.signature());
save(fullfile(sim,'runs/v4/development',sprintf('crossing_t2_w%d_%s.mat',width,label)),'report');
end
