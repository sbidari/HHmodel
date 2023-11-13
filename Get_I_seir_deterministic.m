function [T,I_bar,S_bar,E_bar,R_bar,H_init,P_R,new_inf] = Get_I_seir_deterministic(Time,H_init,P_R,options,Step,dQ,VectN,Vect,nVect,nVectN,NGrid,nTicker,tickGrid,StopProb,E1,E2,DemGrid,alpha,beta,Inf_Ext,g,sigma,Vacc_rate,kB,kV,kL,kR,TB,TV,TL,TR,seed_ext)
nVectS = nVect(1,:);  nVectE = nVect(2,:);  nVectI = nVect(3,:); nVectR = nVect(4,:); 
S_bar = zeros(1,Time/Step); E_bar = zeros(1,Time/Step); I_bar = zeros(1,Time/Step); 
R_bar = zeros(1,Time/Step);     T = zeros(1,(Time/Step)+1);
new_inf = zeros(1,Time/Step);   
t=0;    Struc_Inf = 0;   loop=1; 
Q_demo=Get_Qdemo_seir(P_R,VectN,Vect,Vacc_rate,kB,kV,kL,kR,TB,TV,TL,TR,StopProb);

while (t<Time)
    [Q_int,Q_ext,Q_inf]=Get_Q_seir(Struc_Inf,beta,g,sigma,kB,kV,kL,kR,alpha,t,VectN,Vect,nVect,DemGrid,seed_ext);
    H_inf1 = (Q_inf)'*H_init;
    if mod(t,dQ) == 1
        Q_demo=Get_Qdemo_seir(P_R,VectN,Vect,Vacc_rate,kB,kV,kL,kR,TB,TV,TL,TR,StopProb);
        % Re-calculate P_R
        LeaveStates = find(nTicker==(kB+kL+kV) | nTicker==(kB+kL+kV+kB));   % States at leaving stage
        nER_Leave = nVectR; % Total exposure level at leaving
        P_R = sum( (nER_Leave(LeaveStates)./nVectN(LeaveStates)).*H_init(LeaveStates)')/sum(H_init(LeaveStates));
    end

    [tt,hh]=ode45(@(t,h)Diff(t,h,Q_int,Q_ext,Q_demo),[t t+Step],H_init,options);
    t=tt(end);      H_init=hh(end,:)';    loop=loop+1;    T(loop)=t;

    H_inf2 = (Q_inf)'*H_init;
    new_inf(t) = sum(nVectE.*(H_inf1'+H_inf2')/2);

    % if mod(t,365) == 1
    %     H_init(nVectI==1 & nVectS==2 & nVectE==0 & nVectR==0 & nTicker==4)=H_init(nVectI==1 & nVectS==2 & nVectE==0 & nVectR==0 & nTicker==4)+1;
    % end

    % Calculate the variables that depend on H using aliased H --- P_R, I_T, I_bar 
    % S_bar(loop-1) = round(sum(nVectS.*H_init')); 
    % E_bar(loop-1) = round(sum(nVectE.*H_init'));
    % I_bar(loop-1) = round(sum(nVectI.*H_init'));
    % R_bar(loop-1) = round(sum(nVectR.*H_init'));
    S_bar(loop-1) = sum(nVectS.*H_init');
    E_bar(loop-1) = sum(nVectE.*H_init');
    I_bar(loop-1) = sum(nVectI.*H_init');
    R_bar(loop-1) = sum(nVectR.*H_init');
    % Next loop calculates demographic class-stratified infectious prevalence
    H_T = 0*NGrid;
    Index=cell(1,length(NGrid));
    I_T=zeros(1,length(NGrid));
    for i=1:length(NGrid)
        Index{i}=find(nVectN==NGrid(i)&nTicker==tickGrid(i));
        H_T(i)=sum(H_init(Index{i}));
        I_T(i)=nVectI(Index{i})*H_init(Index{i})/(H_T(i)*NGrid(i));
    end
    Struc_Inf = E1*Inf_Ext*E2*I_T';    
end
T = T(1:end-1);
end

function [dH]=Diff(t,h,Q_int,Q_ext,Q_demo)
dH = (Q_demo+Q_int+Q_ext)'*h;
end
