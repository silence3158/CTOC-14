function pole=v3NominalPole(t,deltaAT,path)
%V3NOMINALPOLE GCRF CIP from bundled arcsecond XYS table (TT), no toolbox.
raw=fileread(path);
epoch=str2double(regexp(raw,'REFEPOCH_JED\s+([\d.]+)','tokens','once'));
step=str2double(regexp(raw,'STEPSIZE\s+([\d.]+)','tokens','once'));
count=str2double(regexp(raw,'NUM_POINTS\s+(\d+)','tokens','once'));
body=extractAfter(raw,'BEGIN TABLE');
values=reshape(sscanf(body,'%f'),3,[]).';
assert(size(values,1)==count,'ctocscreen:v3:xys','XYS table length mismatch.');
% TT = UTC + (TAI-UTC) + 32.184 s. Avoid adding small offsets to a large JD.
jd0=juliandate(datetime(2035,1,1,12,0,0));
u=(jd0-epoch+(t(:)+deltaAT+32.184)/86400)/step;
base=floor(u)-4;
assert(all(base>=0 & base+9<count),'ctocscreen:v3:xys','XYS range exceeded.');
xy=zeros(numel(t),2);
for j=0:9
 weight=ones(size(u));
 for k=0:9
  if k~=j, weight=weight.*((u-base-k)/(j-k)); end
 end
 xy=xy+weight.*values(base+j+1,1:2);
end
xy=xy*(pi/(180*3600));
pole=[xy,sqrt(1-sum(xy.^2,2))];
end
