function [accept,probability,draw]=v3SAAccept(oldCost,newCost,temperature,stream)
%V3SAACCEPT Call only for complete feasible candidates; units km/s.
assert(isfinite(oldCost)&&isfinite(newCost)&&temperature>0);
probability=exp(min(0,(oldCost-newCost)/temperature)); draw=rand(stream); accept=draw<probability;
end
