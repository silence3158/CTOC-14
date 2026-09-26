function key=controlKey(q)
%CONTROLKEY Exact physical controls, independent of the witness plan.
key=sprintf('%.17g,',[q.x0(:);q.T;numel(q.tau);q.tau(:);q.u(:)]);
end
