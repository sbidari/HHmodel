function [E,NGrid,tickGrid,DemGrid] = State_to_Class_Mat_With_Elders(kB,kV,kL,kR,TB,TV,TL,ClassBds,StopProb,points)
% This version adds two parents when children haven't left yet, and two
% elders when they have

% State_to_Class_Mat calculates matrix E where E(T,c) is the expected
% number of people in age class c in a household in demographic state s.
% kB,kL,kB,kR are number of tickers for each transition. ClassBds is a
% list of lower and upper bounds for the classes (with consecutive classes
% touching at boundaries). StopProb is the list of probabilities of not
% having more kids at a certain household size, with first entry always
% zero since we never get here.

maxN = find(StopProb==1,1); % Can't have more kids after this
if isempty(maxN)
    StopProb(end)=1;
    maxN=length(StopProb);
end

% Now calculate set of achievable demographic states
[NGrid,tickGrid]=meshgrid(2:maxN,1:kB+kV+kL+kB+kR);
NGrid = reshape(NGrid,1,numel(NGrid)); 
tickGrid = reshape(tickGrid,1,numel(NGrid));
r=find(NGrid==2 & (tickGrid>kB & tickGrid<=kB+kV+kL+kB)); % Need kids in these states
r=[r find(NGrid==maxN & (tickGrid>kB+kV+kL & tickGrid<=kB+kV+kL+kB))]; % At least one kid has left by now
r=[r find(NGrid==maxN & tickGrid<=kB)]; % Can't have more kids after this
r=[r find(NGrid>2 & tickGrid>kB+kV+kL+kB)]; % Only adults left at this point
NGrid(r) = []; tickGrid(r) = []; % Get rid of bad states
DemGrid = NGrid*(kB+kV+kL+kB+kR)+tickGrid; % List of possible demographic states
KidStateList=find(NGrid>2);

E = zeros(length(DemGrid),length(ClassBds)-1);

for i=1:length(KidStateList) % Only do demo states with kids
    T=KidStateList(i); % Get index of state
    for c=1:length(ClassBds)-2 % Assume kids never hit "elder" class
        [E(T,c),~,~] = State_to_Class(ClassBds(c),ClassBds(c+1),StopProb,NGrid(T),kB,kV,kL,TB,TV,TL,tickGrid(T),30*365,points);
    end
end
E(find(tickGrid<=kB+kV+kL+kB),end-1)=E(find(tickGrid<=kB+kV+kL+kB),end-1)+2; % Add parents % contains both children before leaving home and parents in this class boundary (20-40) years
E(find(tickGrid>kB+kV+kL+kB),end)=2; % Add elders % only old people in this class (older than 40 yrs)

end