function [T,I_bar,S_bar,E_bar,R_bar,H_prev,P_R,new_inf] = Get_I_seir_stochastic(Time,H_prev,P_R,~,Step,dQ,VectN,Vect,nVect,nVectN,NGrid,nTicker,tickGrid,StopProb,E1,E2,DemGrid,alpha,beta,Inf_Ext,g,sigma,Vacc_rate,kB,kV,kL,kR,TB,TV,TL,TR,seed_ext)
nVectS = nVect(1,:);  
nVectE = nVect(2,:);  nVectI = nVect(3,:); nVectR = nVect(4,:);
S_bar = zeros(1,Time/Step); E_bar = zeros(1,Time/Step); I_bar = zeros(1,Time/Step);
R_bar = zeros(1,Time/Step); 
new_inf = zeros(1,Time/Step);   Struc_Inf = 0;
Q_demo=Get_Qdemo_seir(P_R,VectN,Vect,Vacc_rate,kB,kV,kL,kR,TB,TV,TL,TR,StopProb);
for j = 1:Time
    [Q_int,Q_ext,Q_inf]=Get_Q_seir(Struc_Inf,beta,g,sigma,kB,kV,kL,kR,alpha,j,VectN,Vect,nVect,DemGrid,seed_ext);

    if mod(j,dQ) == 1
        Q_demo=Get_Qdemo_seir(P_R,VectN,Vect,Vacc_rate,kB,kV,kL,kR,TB,TV,TL,TR,StopProb);
    end

    Q = (Q_demo+Q_int+Q_ext)';  D=diag(Q);
    P_any = 1 -exp(D*Step);    N_any = binornd(H_prev,P_any);  %Since diagonal contains sum of all rates,
    % N_any = create_brand(H_prev,P_any);
    % this gives all possiblle out events that occur in one step
    Q = Q - diag(diag(Q));  Q_p = Q./sum(Q);   % n*n multinomial probability for individual events
    events = mnrnd(N_any,full(Q_p'));
    change = sum(events)' - N_any;      % compute change for each step using : change = in - out
    H_new = H_prev+ change;

    Qi = Q_p';  Qi(Q_inf<=0)=0;
    dd = ones(length(Qi),1) - sum(Qi,2);    Qi = diag(dd) + Qi;
    events_inf = mnrnd(N_any,full(Qi));
    change_inf = sum(events_inf)' - N_any;
    new_inf(j) = sum(nVectE.*change_inf');

    % Calculate the variables that depend on H using aliased H --- P_R, I_T, I_bar
    S_bar(j) = round(sum(nVectS.*H_new'));
    E_bar(j) = round(sum(nVectE.*H_new'));
    I_bar(j) = round(sum(nVectI.*H_new'));
    R_bar(j) = round(sum(nVectR.*H_new'));
    % Next loop calculates demographic class-stratified infectious prevalence
    H_T = 0*NGrid;
    Index=cell(1,length(NGrid));
    I_T=zeros(1,length(NGrid));
    for i=1:length(NGrid)
        Index{i}=find(nVectN==NGrid(i)&nTicker==tickGrid(i));
        H_T(i)=sum(H_new(Index{i}));
        I_T(i)=nVectI(Index{i})*H_new(Index{i})/(H_T(i)*NGrid(i));
    end
    Struc_Inf = E1*Inf_Ext*E2*I_T';
    % Re-calculate P_R
    LeaveStates = find(nTicker==(kB+kL+kV) | nTicker==(kB+kL+kV+kB));   % States at leaving stage
    nER_Leave = nVectR; % Total exposure level at leaving
    P_R = sum( (nER_Leave(LeaveStates)./nVectN(LeaveStates)).*H_prev(LeaveStates)')/sum(H_prev(LeaveStates));
    H_prev = H_new;
end
T= [1:Time]';
end
