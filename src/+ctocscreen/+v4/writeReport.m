function writeReport(folder,result,trace,eph)
%WRITEREPORT First-round evidence, target table, and physical trajectory plots.
a=result.verification; c=result.manifest.config; s=result.stats;
fid=fopen(fullfile(folder,'report.md'),'w','n','UTF-8'); assert(fid>0); clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# V4 cold-start experiment\n\n');
fprintf(fid,'Seed: %d. Search budget: %.1f s. Measured search: %.3f s.\n\n',c.seed,c.budget_s,s.search_seconds);
fprintf(fid,'Independent target visits: **%d/35**. Raw total delta-V: **%.12f km/s**. Complete physical pass: **%d**.\n\n',a.visit_count,a.total_dv_km_s,a.passed);
fprintf(fid,'Height lower bound: %.9f km. Initial orbit pass: %d. Duration: %.3f s.\n\n', ...
 a.min_altitude_lower_km,a.initial_passed,result.best.q.T);
fprintf(fid,'Source signature unchanged: %d. Empty historical input: 1. Official alignment verified: 0.\n\n',s.source_unchanged);
fprintf(fid,'Full-history B calls: %d. Shared-arc calls: %d. Structure calls: %d. Joint iterations: %d.\n\n', ...
 s.full_calls,s.shared_calls,s.structure_calls,s.joint_iterations);
fprintf(fid,'Search pool threshold: 6.1 km/s. User notification threshold: <8 km/s with independent 35/35.\n\n');
fprintf(fid,'| Target | Independent distance (km) | Witness time (s) | Visited |\n|---|---:|---:|---:|\n');
for j=1:35, fprintf(fid,'| %d | %.12g | %.9f | %d |\n',j,a.distance_km(j),a.witness_times_s(j),a.distance_km(j)<=1); end
fprintf(fid,'\n## Controls\n\n| Pulse | Time (s) | dv-x | dv-y | dv-z | Magnitude (km/s) |\n|---|---:|---:|---:|---:|---:|\n');
q=result.best.q;
for k=1:numel(q.tau), fprintf(fid,'| %d | %.9f | %.12g | %.12g | %.12g | %.12g |\n',k,q.tau(k),q.u(k,:),norm(q.u(k,:))); end
fprintf(fid,'\nThe report distinguishes actual physical results from solver residuals. Missing visits make this an incomplete mission result.\n');
samples=zeros(0,7);
for k=1:numel(trace.arcs)
 tt=linspace(trace.arc_start(k),trace.arc_end(k),max(20,ceil((trace.arc_end(k)-trace.arc_start(k))/300)));
 yy=deval(trace.arcs{k},tt); samples=[samples;tt.',yy(1:6,:).']; %#ok<AGROW>
end
save(fullfile(folder,'trajectory_samples.mat'),'samples');
try
 fig=figure('Visible','off','Color','w');
 if ~isempty(samples), plot3(samples(:,2),samples(:,3),samples(:,4),'LineWidth',1); hold on; end
 rt=eph.states0(:,1:3); scatter3(rt(:,1),rt(:,2),rt(:,3),18,'filled');
 if ~isempty(samples), legend('Inspector trajectory','Targets at initial epoch','Location','best'); end
 axis equal; grid on; xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
 title(sprintf('V4: %d/35, raw delta-V %.3f km/s',a.visit_count,a.total_dv_km_s));
 exportgraphics(fig,fullfile(folder,'trajectory.png'),'Resolution',140); close(fig);
 fig=figure('Visible','off','Color','w');
 if ~isempty(q.tau), stairs([0;q.tau;q.T],[0;cumsum(vecnorm(q.u,2,2));a.total_dv_km_s],'LineWidth',1.3); end
 grid on; xlabel('Mission time (s)'); ylabel('Cumulative raw delta-V (km/s)');
 exportgraphics(fig,fullfile(folder,'delta_v.png'),'Resolution',140); close(fig);
catch err
 fprintf('V4 plot generation: %s\n',err.message);
end
end
