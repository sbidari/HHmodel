function [Q_int,Q_ext,Q_inf]=Get_Q_seir(Struc_Inf,beta_0,g,sigma,kB,kV,kL,kR,alpha,t,VectN,Vect,nVect,DemGrid,seed)

VectS=Vect(1,:);    VectE=Vect(2,:);    VectI=Vect(3,:);    VectR=Vect(4,:);
nVectS=nVect(1,:);  nVectE=nVect(2,:);  nVectI=nVect(3,:);  nVectR=nVect(4,:);

% Add a seasonal forcing term for beta 
beta = beta_0*(1 + alpha * cos((2*pi*t)/365) );
% SIR classes without external transmission.

Q_int_short=sparse(1,1,0,length(VectN),length(VectN)); % Transition matrix for internal epidemic
Q_inf_short=sparse(1,1,0,length(VectN),length(VectN)); % Transition matrix for counting incidence 
for LeaveStates=1:length(VectN)
    R=(beta*VectI(LeaveStates)/(VectN(LeaveStates)-1))*VectS(LeaveStates) + sigma*VectE(LeaveStates) + g*VectI(LeaveStates) ; 
    Q_int_short(LeaveStates,LeaveStates)=-R;            
    
    R=(beta*VectI(LeaveStates)/(VectN(LeaveStates)-1))*VectS(LeaveStates) ; % Infection
    Q_inf_short(LeaveStates,LeaveStates)=-R;
    if R>0
        q=find(VectS==VectS(LeaveStates)-1 & VectE==VectE(LeaveStates)+1 & VectI==VectI(LeaveStates) & VectR==VectR(LeaveStates));
        Q_int_short(LeaveStates,q)=Q_int_short(LeaveStates,q)+R; % Add to transition matrix
        Q_inf_short(LeaveStates,q)=Q_inf_short(LeaveStates,q)+R;
    end

    R=sigma*VectE(LeaveStates); % Infectious
    if R>0
        q=find(VectS==VectS(LeaveStates) & VectE==VectE(LeaveStates)-1 & VectI==VectI(LeaveStates)+1 & VectR==VectR(LeaveStates));
        Q_int_short(LeaveStates,q)=Q_int_short(LeaveStates,q)+R;
    end

    R=g*VectI(LeaveStates); % Recovery
    if R>0
        q=find(VectS==VectS(LeaveStates) & VectE==VectE(LeaveStates) & VectI==VectI(LeaveStates)-1 & VectR==VectR(LeaveStates)+1);
        Q_int_short(LeaveStates,q)=Q_int_short(LeaveStates,q)+R;
    end
end

Tickers=(kB)+(kV)+(kL)+(kB)+(kR); % Total number of possible tick values
L=size(Vect,2); 

% First Increase the Basic Sizes
% L is size of standard system, j is Ticker position.
% So nVectX([1:L]) is the original system, nVectX(L+[1:L]) is the 2nd
% ticker etc etc.
for i=Tickers:-1:1 % Puts in full range of ticker classes for each epidemic state
    j=i-1;
    nVectN(j*L+[1:L])=VectN;
    nTicker(j*L+[1:L])=i;
    Q_int(j*L+[1:L],j*L+[1:L])=Q_int_short; % Add a given ticker state to transmission rate matrix
    Q_inf(j*L+[1:L],j*L+[1:L])=Q_inf_short;
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

nVectN(r)=[];  nTicker(r)=[];  Q_int(r,:)=[];  Q_int(:,r)=[]; Q_inf(r,:)=[];  Q_inf(:,r)=[];

%% We now calculate eM, the transition matrix of external infection events.

Struc_Inf = Struc_Inf*(1 + alpha * cos((2*pi*t)/365) );
T = nVectN*(kB+kV+kL+kB+kR)+nTicker; % Full demo state of each demo-epi state

if size(Struc_Inf)==1
    Struc_Inf=Struc_Inf*ones(1,length(DemGrid));
end

Q_ext=sparse(1,1,0,length(nVectN),length(nVectN)); % eM is external infection
for LeaveStates=1:length(nVectN)  
    
    R=Struc_Inf(find(DemGrid==T(LeaveStates))) + seed;
    
    R=R*nVectS(LeaveStates)/nVectN(LeaveStates); %%% This assumes external is per person %%%
    
    if R>0
        Q_ext(LeaveStates,LeaveStates)=-R;  % This doesnt include recovery since recovery is already included in Q_int
        Q_inf(LeaveStates,LeaveStates)=Q_inf(LeaveStates,LeaveStates)-R;
        
        q=find(nVectS==nVectS(LeaveStates)-1 & nVectE==nVectE(LeaveStates)+1 & nVectI==nVectI(LeaveStates) & nVectR==nVectR(LeaveStates) & nTicker==nTicker(LeaveStates));
        Q_ext(LeaveStates,q)=Q_ext(LeaveStates,q)+R;
        Q_inf(LeaveStates,q)=Q_inf(LeaveStates,q)+R;
    end
end

end 




