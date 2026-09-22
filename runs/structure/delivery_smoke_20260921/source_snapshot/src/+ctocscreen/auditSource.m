function audit=auditSource(problem,scenarioPath,gravityPath)
%AUDITSOURCE Offline read-only comparison, never starts ATK.
doc=xmlread(scenarioPath); satellites=doc.getElementsByTagName('Satellite');
states=nan(35,6); epochs=strings(35,1);
tags={'PositionX','PositionY','PositionZ','VelocityX','VelocityY','VelocityZ'};
for k=0:satellites.getLength-1
 node=satellites.item(k); name=char(node.getAttribute('Name'));
 id=sscanf(name,'Target%d');
 if isempty(id)||id<1||id>35, continue; end
 for j=1:6
  item=node.getElementsByTagName(tags{j}).item(0);
  states(id,j)=str2double(char(item.getTextContent))/1000;
 end
 epochs(id)=string(node.getElementsByTagName('OrbEpoch').item(0).getTextContent);
end
assert(all(isfinite(states),'all'),'Missing original target states.');
assert(all(epochs==problem.epoch_utc),'Original target epoch mismatch.');
err=abs(states-problem.states0);
assert(max(err(:,1:3),[],'all')<=5.1e-7,'CSV position extraction mismatch.');
assert(max(err(:,4:6),[],'all')<=5.1e-10,'CSV velocity extraction mismatch.');
fid=fopen(gravityPath); cleanup=onCleanup(@()fclose(fid)); values=fscanf(fid,'%f',4);
assert(abs(values(3)/1e9-problem.mu_km3_s2)<1e-8 && ...
 abs(values(4)/1000-problem.re_km)<1e-8,'Gravity constants mismatch.');
audit=struct('status','passed','max_position_component_error_km',max(err(:,1:3),[],'all'), ...
 'max_velocity_component_error_km_s',max(err(:,4:6),[],'all'), ...
 'frame','ATK_native_inertial','frame_note','Native Cartesian axes retained; ICRF/J2000 identification not certified.', ...
 'scenario_path',scenarioPath,'gravity_path',gravityPath);
end
