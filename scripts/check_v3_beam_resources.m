function check_v3_beam_resources()
%CHECK_V3_BEAM_RESOURCES Preserve scarce mission time in partial beam labels.
sim=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(sim,'src'));
c=ctocscreen.v3Defaults();
base=struct('state',[7000 0 0 0 7.5 0],'visited',[true(10,1);false(25,1)], ...
 't',600000,'J',3,'estimate',3,'last_gain',1);
late=base; cheaper=base; cheaper.J=2.9; cheaper.estimate=2.9; cheaper.t=590000;
early=base; early.t=200000; early.J=3.2; early.estimate=3.2;
coast=base; coast.last_gain=0; coast.t=100000; coast.state(2)=3000;
selected=ctocscreen.v3SelectBeam({late,cheaper,early,coast},3,c);
assert(any(cellfun(@(n)n.t==early.t,selected)));
assert(any(cellfun(@(n)n.J==cheaper.J,selected)));
assert(any(cellfun(@(n)n.last_gain==0,selected)));
% An already complete mission has no unfinished time resource to preserve.
late.visited(:)=true; early.visited(:)=true;
selected=ctocscreen.v3SelectBeam({late,early},1,c); assert(selected{1}.J==late.J);
fprintf('BEAM_RESOURCES_CHECK: early prefix, cheap prefix, zero-visit opportunity retained.\n');
end
