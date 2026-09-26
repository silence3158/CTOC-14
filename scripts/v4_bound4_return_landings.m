function [path,returns,stopped,legs,hits]=v4_bound4_return_landings(path,burns,returns,stopped,depth,limit)
%V4_BOUND4_RETURN_LANDINGS Count actual landings only; cascade at the limit.
% burns is aligned with path. returns/stopped are indexed by global node ID.
% legs columns: from, to, removed burns, destination return count.
legs=zeros(0,4); hits=[];
while numel(path)>1
 M=burns(end);
 if M<depth, pos=1; else, pos=find(burns(1:end-1)<=M-depth,1,'last'); end
 assert(~isempty(pos),'No four-impulse ancestor.');
 removed=M-burns(pos);
 assert(removed==min(depth,M),'Incorrect rollback burn count.');
 destination=path(pos);
 returns(destination)=returns(destination)+1;
 legs(end+1,:)=[path(end),destination,removed,returns(destination)]; %#ok<AGROW>
 path=path(1:pos); burns=burns(1:pos);
 if returns(destination)<limit, return; end
 assert(~stopped(destination),'Landed on an already stopped node.');
 stopped(destination)=true; hits(end+1)=destination; %#ok<AGROW>
 if isscalar(path), path=[]; return; end
end
end
