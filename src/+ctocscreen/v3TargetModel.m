function [model,audit,p] = v3TargetModel()
%V3TARGETMODEL Offline audit and explicit nominal Earth orientation model.
% Native XML axes are preserved numerically. GCRF interpretation is explicit
% and provisional until competition-engine alignment has been established.
p=ctocscreen.loadProblem();
g=p.source_audit.gravity_path;
raw=sscanf(fileread(g),'%f'); rows=reshape(raw(5:end),4,[]).';
c20=rows(rows(:,1)==2 & rows(:,2)==0,3);
assert(isscalar(c20) && c20<0,'ctocscreen:v3:gravity','Missing C20.');
model=struct('dynamics_id','central_j2','mu',raw(3)/1e9,'re',raw(4)/1000, ...
 'c20_normalized',c20,'j2',-sqrt(5)*c20,'horizon_s',864000, ...
 'epoch_utc','2035-01-01 12:00:00','frame','ATK_native_assumed_GCRF', ...
 'orientation','Bundled IAU2006 XYS, order-9 Lagrange in TT; zero polar motion/dCIP; TAI-UTC held at 37 s', ...
 'official_alignment_verified',false);
atk=fileparts(fileparts(fileparts(g)));
leap=fileread(fullfile(atk,'AstroData','LeapSecond.dat'));
tokens=regexp(leap,'(?m)^\s*(\d{5})\s+([\d.]+)\s+\d{4}','tokens');
assert(~isempty(tokens),'ctocscreen:v3:leap','No leap-second data.');
model.tai_minus_utc_s=str2double(tokens{end}{2});
assert(model.tai_minus_utc_s==37,'ctocscreen:v3:leap','Review new leap data.');
xysPath=fullfile(atk,'AstroData','IERS-conventions','2010','IAU2006_XYS.dat');
model.pole_times=(0:300:model.horizon_s).';
model.pole=ctocscreen.v3NominalPole(model.pole_times,model.tai_minus_utc_s,xysPath);
mid=model.pole_times(1:end-1)+150;
direct=ctocscreen.v3NominalPole(mid,model.tai_minus_utc_s,xysPath);
approx=interp1(model.pole_times,model.pole,mid,'linear');
approx=approx./vecnorm(approx,2,2);
model.pole_midpoint_error_rad=max(vecnorm(direct-approx,2,2));
assert(model.pole_midpoint_error_rad<1e-8,'ctocscreen:v3:poleAccuracy','Pole interpolation inaccurate.');
audit=p.source_audit;
audit.c20_normalization='fully normalized; J2=-sqrt(5)*C20';
audit.status='input_values_checked_frame_alignment_pending';
audit.frame_interpretation=model.frame;
audit.frame_evidence='Bundled ATK frame manual describes GCRF; XML contains no explicit native-axis label. No rotation applied to CSV.';
audit.orientation_assumptions=model.orientation;
audit.official_alignment_verified=false;
eop=fileread(fullfile(atk,'AstroData','EOP-All.txt'));
dates=regexp(eop,'(?m)^\s*(\d{4})\s+(\d{1,2})\s+(\d{1,2})\s+\d{5}','tokens');
audit.eop_last_date=strjoin(dates{end},'-');
audit.future_eop_note='2035 is outside supplied EOP coverage. Zero-EOP nominal model is NOT an asserted ATK fallback.';
doc=xmlread(p.source_audit.scenario_path); sats=doc.getElementsByTagName('Satellite');
ids=[];
for k=0:sats.getLength-1
 n=sats.item(k); id=sscanf(char(n.getAttribute('Name')),'Target%d');
 if isempty(id) || id<1 || id>35, continue; end
 ids(end+1)=id; %#ok<AGROW>
 orbit=n.getElementsByTagName('Orbit').item(0);
 assert(strcmp(char(n.getElementsByTagName('CentralBody').item(0).getTextContent),'Earth'), ...
  'ctocscreen:v3:scenario','Target central body must be Earth.');
 keys={'OrbPropType','GravityModel','MaxDegree','MaxOrder','UseDrag','UseSRP','UseRelativity'};
 expected=[2 0 2 0 0 0 0];
 for j=1:numel(keys)
  value=str2double(char(orbit.getElementsByTagName(keys{j}).item(0).getTextContent));
  assert(value==expected(j),'ctocscreen:v3:scenario','Unexpected target force settings.');
 end
 bodies=orbit.getElementsByTagName('ThirdBodyGrav');
 for j=0:bodies.getLength-1
  body=bodies.item(j); code=str2double(char(body.getAttribute('CbCode')));
  enabled=str2double(char(body.getElementsByTagName('UseGravity').item(0).getTextContent));
  assert(code==2 || enabled==0,'ctocscreen:v3:scenario','Unexpected enabled third body.');
 end
end
assert(isequal(sort(ids),1:35),'ctocscreen:v3:scenario','Missing or duplicate original targets.');
audit.target_force_settings_checked=true;
audit.user_frame_statement='2026-09-23: extracted directly from XML; axes not recalled; likely Earth-centered. This does not identify axes.';
end
