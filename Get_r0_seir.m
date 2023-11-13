function [r] = Get_r0_seir(Inf_Ext,beta,g,sigma,popH,kB,kL,kV,kR,TB,TL,TV,TD,Distrib_Children,StopProb,Vacc_rate,DemGrid,NGrid,tickGrid,E1,m_run,seed)

% Calculate disease-free equilibrium distribution and internal transmission matrix
[~,H_Eq,~,VectN,Vect,nVect,nTicker,nVectN,E1,E2] = Get_Eq_Demography_seir(popH,TB,kB,TL,kL,kV,TV,TD,kR,Distrib_Children,StopProb,Vacc_rate,NGrid,tickGrid,E1);

%set alpha amplitude of sasonality as 0 for now
H_Eq = H_Eq./sum(H_Eq);
[Q_int,~,~] = Get_Q_seir(0,beta,g,sigma,kB,kV,kL,kR,0,0,VectN,Vect,nVect,DemGrid,seed); 

% Compute H_T
H_T = 0*NGrid;
Index=cell(1,length(NGrid));
for i=1:length(NGrid)
    Index{i}=find(nVectN==NGrid(i)&nTicker==tickGrid(i));
    H_T(i)=sum(H_Eq(Index{i}));
end

T = nVectN*(kB+kL+kV+kB+kR)+nTicker; % Indices of full demographic states
F_D=H_Eq; F_D(F_D<0)=0; F_D(nVect(2,:)==0 & nVect(3,:)==1 & nVect(4,:)==0)=1e-5;  % add a bit of infection
Q_int_long=Q_int;
m=find(nVect(2,:)+nVect(3,:)+nVect(4,:)==0);        %RECONSIDER THIS!  ok!
for i=1:length(m)
    Q_int_long(m(i),m(i))=0; % Get rid of outflow so that we can study growth in infected states
end

% In this section we calculate exponential growth rate r=\dot{I}/I.
% We repeatedly run the linearised equations forward over short distances,
% scaling down the additional infection at the end of each run. In this way
% we get the quasi-equilibrium distribution of infection in the early
% stages of the epidemic, which is the eigenvector associated with the
% eigenvalue r, which is the growth rate of the system.
flag=1;
attempts=0;
while flag
    [t, X] = ode45( @(t,p) ODEs(t,p,nVect, nTicker, Q_int_long, DemGrid, T, H_T, E1, E2, Inf_Ext,kB+kL+kV+kB+kR),...
        [0 m_run], F_D, odeset('RelTol',1e-8, 'AbsTol', 1e-10,'NonNegative',[1:length(F_D)]) );
    FFD=X(end,:);

    FFD(nVect(2,:)+nVect(3,:)+nVect(4,:)>0) = 1e-5*FFD(nVect(2,:)+nVect(3,:)+nVect(4,:)>0)/sum(FFD(nVect(2,:)+nVect(3,:)+nVect(4,:)>0));
    dX=ODEs(0,X(end,:)',nVect, nTicker, Q_int_long, DemGrid, T, H_T, E1, E2, Inf_Ext,kB+kL+kV+kB+kR);
    Xnow=X(end,:)';
    R=dX./Xnow; R=R((Xnow'>0)&(nVect(2,:)+nVect(3,:)>0));   %INLCUDE BOTH E&I
    attempts=attempts+1;
    sXn=sum(Xnow(nVect(2,:)+nVect(3,:)>0)); sFD=sum(F_D(nVect(2,:)+nVect(3,:)>0));
    fprintf(1,'Attempt %d) X=%g R=%g var(R)=%g [sum(end)=%g sum(start)=%g %g]\n',attempts,mean(Xnow),mean(R),var(R),sXn,sFD,log(sXn/sFD)/t(end));
    if var(R)<1e-7||attempts>99
        flag=0;
        if attempts>99
            R=log(sXn/sFD)/t(end);
        end
    end
    if sXn<sFD*5
        FFD(nVect(2,:)+nVect(3,:)+nVect(4,:)>0)=1e-5*FFD(nVect(2,:)+nVect(3,:)+nVect(4,:)>0)/sum(FFD(nVect(2,:)+nVect(3,:)>0));
        if sXn<sFD*0.1 && attempts>1 && flag
            flag=0;
            R=log(sXn/sFD)/t(end);
        end
        if sXn<sFD*10 && attempts>10 && flag
            flag=0;
            R=log(sXn/sFD)/t(end);
        end
    end
    F_D=FFD;
end
r=mean(R); % Converged per capita growth rate

end

function [dX] = ODEs(t, Z, nVect, nTicker, Q_int, DemGrid, T, H_T, E1, E2, D,sumK)
% This does the master equations except with diagonal set to zero, so
% probability mass accumulates in states rather than flowing out. This
% allows us to measure the growth in infection in the early stages of an
% outbreak.

nVectS=nVect(1,:);  nVectE=nVect(2,:);  nVectI=nVect(3,:);   nVectR=nVect(4,:);
nVectN=sum(nVect,1); 
I_T=zeros(1,length(H_T));

for i=1:length(nVectN)
    State_i = find(DemGrid==T(i));
    I_T(State_i) = I_T(State_i) + (nVectI(i)/nVectN(i))*Z(i);
end
I_T=I_T./H_T;
Struc_Inf = E1*D*E2*I_T'; % Age-structured ext. inf.

Q_ext=sparse(1,1,0,length(nVectN),length(nVectN));
for n=1:length(nVectN) % Builds matrix of external infection rates
    if (nVect(2,n)+nVect(3,n)+nVect(4,n))==0  %% as assuming infection is rare only to all S's (branching assumption - Joe)
        R=Struc_Inf(find(DemGrid==nVectN(n)*(sumK)+nTicker(n)))/nVectN(n);
        R=R*nVectS(n); %%% This assumes external is per person %%%
        if R>0
            q=find(nVectS==nVectS(n)-1 & nVectE==nVectE(n)+1 & nVectI==nVectI(n) & nVectR==nVectR(n) & nTicker==nTicker(n));
            Q_ext(n,q)=Q_ext(n,q)+R;
        end
    end
end
dX=(Q_int+Q_ext)'*Z;
end
