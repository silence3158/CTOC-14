function [C,S] = stumpff(z)
%STUMPFF Stable universal-variable functions for a real scalar.
if abs(z)<1e-4
 C=0.5-z/24+z^2/720-z^3/40320+z^4/3628800;
 S=1/6-z/120+z^2/5040-z^3/362880+z^4/39916800;
elseif z>0
 w=sqrt(z); C=2*(sin(w/2)/w)^2; S=(w-sin(w))/w^3;
else
 w=sqrt(-z); C=2*(sinh(w/2)/w)^2; S=(sinh(w)-w)/w^3;
end
end