function check_v4_return_landings
% Synthetic scheduler checks only; no orbital results or historical inputs.
R=zeros(1,12); S=false(1,12);
% Five impulses beyond B (node 5): land at node 6, leave B unchanged.
[p,r,s,legs,h]=v4_bound4_return_landings(1:10,0:9,R,S,4,30);
assert(isequal(p,1:6)&&r(5)==0&&r(6)==1&&~any(s)&&isempty(h)&&legs(1,3)==4);
% The 29th actual return stays; the 30th cascades another four impulses.
R(5)=28;
[p,r,~,~,h]=v4_bound4_return_landings(1:9,0:8,R,S,4,30);
assert(isequal(p,1:5)&&r(5)==29&&isempty(h));
[p,r,s,legs,h]=v4_bound4_return_landings(1:9,0:8,r,S,4,30);
assert(isequal(p,1)&&r(5)==30&&r(1)==1&&s(5)&&isequal(h,5)&&size(legs,1)==2);
% Both landing points hit the threshold; retire root.
R=zeros(1,12); R([1 5])=29;
[p,r,s,legs,h]=v4_bound4_return_landings(1:9,0:8,R,S,4,30);
assert(isempty(p)&&all(r([1 5])==30)&&all(s([1 5]))&&isequal(h,[5 1])&&all(legs(:,3)==4));
% A shallow branch removes only available impulses and increments root once.
[p,r,~,legs,~]=v4_bound4_return_landings(1:3,0:2,zeros(1,12),S,4,30);
assert(isequal(p,1)&&r(1)==1&&sum(r)==1&&legs(1,3)==2);
% Zero-impulse coasting nodes do not consume the four-impulse depth.
[p,r,~,legs,~]=v4_bound4_return_landings(1:7,[0 0 1 2 2 3 4],zeros(1,12),S,4,30);
assert(isequal(p,1:2)&&r(2)==1&&legs(1,3)==4);
% Disabled threshold retains legacy single-jump behavior.
[p,r,s,legs,h]=v4_bound4_return_landings(1:9,0:8,R,S,4,Inf);
assert(isequal(p,1:5)&&r(5)==30&&~any(s)&&isempty(h)&&size(legs,1)==1);
fprintf('PASS: six synthetic actual-landing scheduler scenarios.\n');
end
