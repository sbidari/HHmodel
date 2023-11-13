function [TR,DiseaseFree,P_R,VectN,Vect,nVect,nTicker,nVectN,E1,E2] = Get_Eq_Demography_seir(popH,TB,kB,TL,kL,kV,TV,TD,kR,Distrib_Children,StopProb,Vacc_rate,NGrid,tickGrid,E1)

% Define demographic parameters
Exp_Children = (0:length(Distrib_Children)-1)*Distrib_Children'; % Expected number of children
TR = TD-TL-TV-TL-TV-TB*Exp_Children; % Reset rate = age of death - expected time with kids
maxN = find(StopProb==1,1); % Can't have more kids after this
if isempty(maxN)
    StopProb(end)=1;
    maxN=length(StopProb);
end

[S,E,I,R]=ndgrid([0:maxN],[0:maxN],[0:maxN],[0:maxN]);
m=find(S+E+I+R<=maxN & S+E+I+R>=2)'; % Impose household size constraintes
VectN=S(m)+E(m)+I(m)+R(m); VectS=S(m); VectE=E(m); VectI=I(m); VectR=R(m);
Vect=[VectS; VectE; VectI; VectR]; % Ordering of all states

Tickers=(kB)+(kV)+(kL)+(kB)+(kR); % Total number of possible tick values
L=size(Vect,2);

for i=Tickers:-1:1 % Puts in full range of ticker classes for each epidemic state
    j=i-1;
    nVect(:,j*L+[1:L])=Vect;
    nVectN(j*L+[1:L])=VectN;
    nTicker(j*L+[1:L])=i;
end

% NOW NEED TO REMOVE SURPLUS POINTS
% if waiting for kids to leave and only size 2
r=find(nVectN==2 & (nTicker>kB & nTicker<=kB+kV+kL+kB));
% if waiting for second kids to leave and size=max
r=[r find(nVectN==max(nVectN) & (nTicker>kB+kV+kL & nTicker<=kB+kV+kL+kB))];
% if waiting for birth and size = max
r=[r find(nVectN==max(nVectN) & nTicker<=kB)];
% if waiting for reset and bigger than size 2
r=[r find(nVectN>2 & nTicker>kB+kV+kL+kB)];
% This syntax just deletes these entries (changing size of array)
nVect(:,r)=[]; nVectN(r)=[];    nTicker(r)=[];

%First start by computing the demographic equilibria
nVectS = nVect(1,:);    nVectI = nVect(3,:);    nVectR = nVect(4,:);


diff = 1;   P_R_new=0;  i=1;
while diff(i) >1e-6
    P_R_old = P_R_new;
    Q_demo=Get_Qdemo_seir(P_R_old,VectN,Vect,Vacc_rate,kB,kV,kL,kR,TB,TV,TL,TR,StopProb);
    [H,~] = eigs(Q_demo',1,'lr');
    H = H/sum(H);
    H(H<1e-12) = 0;
    H_Eq = round(H*popH);
    % Calculate P_R
    LeaveStates=find(nTicker==(kB+kL+kV) | nTicker==(kB+kL+kV+kB)); % States at leaving stage
    nER_Leave = nVectR; % Total exposure level at leaving
    P_R_new = sum( (nER_Leave(LeaveStates)./nVectN(LeaveStates)).*H_Eq(LeaveStates)')/sum(H_Eq(LeaveStates));
    i = i+1;
    diff(i) = abs(P_R_new - P_R_old);
end
P_R = P_R_new;

DiseaseFree=round(H_Eq);
% Next loop calculates demographic class-stratified infectious prevalence
H_T = 0*NGrid;
Index=cell(1,length(NGrid));
for i=1:length(NGrid)
    Index{i}=find(nVectN==NGrid(i)&nTicker==tickGrid(i));
    H_T(i)=sum(H_Eq(Index{i}));
end

% The following loop makes sure that all the rows of E sum to the household
% size. This is necessary because we calculated it using Monte Carlo
% integration, which will not be completely accurate.
for i=1:length(NGrid)
    if NGrid(i)>2
        E1(i,end-1)=E1(i,end-1)-2;
        E1(i,:)=(NGrid(i)-2)*E1(i,:)/sum(E1(i,:));
        E1(i,end-1)=E1(i,end-1)+2;
    end
end
PQ_Sum = H_T*E1;
E2=H_T.*E1'./PQ_Sum';

end

