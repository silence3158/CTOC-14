classdef TestV3Preprocessing < matlab.unittest.TestCase
 % Set CTOC_V3_TEST_CACHE to a completed cache filename before running.
 properties
  Cache
  Filename
 end
 methods(TestClassSetup)
  function loadCache(test)
   test.Filename=getenv('CTOC_V3_TEST_CACHE');
   test.assertTrue(isfile(test.Filename),'Set CTOC_V3_TEST_CACHE first.');
   test.Cache=ctocscreen.v3LoadTargetEphemeris(test.Filename);
  end
 end
 methods(Test)
  function initialStates(test)
   [r,v]=ctocscreen.v3QueryTargets(test.Cache,[],0);
   test.verifyEqual([r v],test.Cache.states0,'AbsTol',1e-11);
  end
  function fullHorizon(test)
   [r,v,a]=ctocscreen.v3QueryTargets(test.Cache,[],[0 864000],'grid');
   test.verifySize(r,[35 3 2]);
   test.verifyTrue(all(isfinite([r(:);v(:);a(:)])));
  end
  function gridAndPairs(test)
   [r,v]=ctocscreen.v3QueryTargets(test.Cache,[35 2 35],[864000 0 12345],'grid');
   [rp,vp]=ctocscreen.v3QueryTargets(test.Cache,[35 2 35],[864000 0 12345]);
   test.verifyEqual([r(1,:,1);r(2,:,2);r(3,:,3)],rp);
   test.verifyEqual([v(1,:,1);v(2,:,2);v(3,:,3)],vp);
  end
  function velocityDerivative(test)
   t=200123; h=0.01;
   [~,v]=ctocscreen.v3QueryTargets(test.Cache,[],t);
   plus=ctocscreen.v3QueryTargets(test.Cache,[],t+h);
   minus=ctocscreen.v3QueryTargets(test.Cache,[],t-h);
   test.verifyLessThan(max(vecnorm((plus-minus)/(2*h)-v,2,2)),1e-7);
  end
  function knotContinuity(test)
   h=test.Cache.targets{1}.step_s; t=100*h; epsT=1e-4;
   [r,v]=ctocscreen.v3QueryTargets(test.Cache,1,[t-epsT t t+epsT]);
   test.verifyLessThan(norm(r(3,:)-r(1,:)-2*epsT*v(2,:)),1e-8);
   test.verifyLessThan(norm(v(3,:)-v(1,:)),1e-5);
  end
  function rejectOutOfRange(test)
   test.verifyError(@()ctocscreen.v3QueryTargets(test.Cache,1,-1),'ctocscreen:v3:timeRange');
   test.verifyError(@()ctocscreen.v3QueryTargets(test.Cache,1,864001),'ctocscreen:v3:timeRange');
  end
  function rejectInvalidInputs(test)
   test.verifyNotEmpty(errorId(@()ctocscreen.v3QueryTargets(test.Cache,36,0)));
   test.verifyNotEmpty(errorId(@()ctocscreen.v3QueryTargets(test.Cache,1,NaN)));
   test.verifyError(@()ctocscreen.v3QueryTargets(test.Cache,[1 2],[0 1 2]),'ctocscreen:v3:shape');
   test.verifyError(@()ctocscreen.v3QueryTargets(test.Cache,1,0,'bad'),'ctocscreen:v3:mode');
  end
  function nominalIsNotOfficial(test)
   test.verifyFalse(test.Cache.model.official_alignment_verified);
   test.verifyError(@()ctocscreen.v3LoadTargetEphemeris(test.Filename,true),'ctocscreen:v3:alignmentPending');
  end
  function alteredSignatureRejected(test)
   cache=test.Cache; cache.signature.sha256{1}='changed';
   file=[tempname '.mat']; cleanup=onCleanup(@()removeCache(file));
   save(file,'cache');
   fid=fopen([file '.sha256'],'w'); fprintf(fid,'%s',ctocscreen.v3FileHash(file)); fclose(fid);
   test.verifyError(@()ctocscreen.v3LoadTargetEphemeris(file),'ctocscreen:v3:staleCache');
  end
  function corruptedCacheRejected(test)
   file=[tempname '.mat']; cleanup=onCleanup(@()removeCache(file));
   fid=fopen(file,'w'); fprintf(fid,'corrupt'); fclose(fid);
   fid=fopen([file '.sha256'],'w'); fprintf(fid,'wrong'); fclose(fid);
   test.verifyError(@()ctocscreen.v3LoadTargetEphemeris(file),'ctocscreen:v3:checksum');
  end
  function equatorialAndPolarForce(test)
   m=test.Cache.model; m.pole(:,:)=repmat([0 0 1],size(m.pole,1),1);
   R=7000; f=1.5*m.j2*(m.re/R)^2;
   a=ctocscreen.v3J2Acceleration(0,[R 0 0;0 0 R],m);
   expected=-m.mu/R^2*[1+f 0 0;0 0 1-2*f];
   test.verifyEqual(a,expected,'AbsTol',1e-16);
  end
  function numericalBudget(test)
   test.verifyTrue(test.Cache.validation.passed);
   test.verifyLessThanOrEqual(test.Cache.validation.max_position_error_m,0.1);
   test.verifyLessThanOrEqual(test.Cache.validation.max_velocity_error_m_s,1e-4);
   test.verifySize(test.Cache.validation.table,[35 8]);
  end
 end
end

function removeCache(file)
if isfile(file), delete(file); end
if isfile([file '.sha256']), delete([file '.sha256']); end
end

function id=errorId(f)
id='';
try
 f();
catch err
 id=err.identifier;
end
end
