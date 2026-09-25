function out=v3Verify(s,eph,c)
%V3VERIFY Export gate; tighten and repeat if a successful witness is near 1 km.
out=ctocscreen.v3Replay(s,eph,c,true); out.stability_repeated=false;
if out.passed && (max(out.distance_km)>0.95 || out.sampled_min_altitude_km<200.1)
 fine=c; fine.max_step_s=min(30,c.max_step_s/2); fine.verify_reltol=max(2.3e-14,c.verify_reltol/1.2);
 second=ctocscreen.v3Replay(s,eph,fine,true);
 out.stability_repeated=true; out.repeat=second;
 out.passed=second.passed;
 if ~second.passed, out.status='tightened_verification_failed'; out.validation_level='proposal'; end
end
end
