function[t,I_bar,S_bar,R_bar,E_bar,daily_new_inf,yearly_infection]=run_epi_1968_2019
% Generates annual measles incidence from 1968 - 2019 assuming changing household
% size according to the household size data and stores the outcome in a mat-file
% The simulation can be run faster using mex file (written in C++) to
% generate random numbers. This is disabled here to prevent issues with
% computer architecture and compilation


% mex create_brand.cpp
% params fixed
popH = 5e4;
seed = 1e-6;
g = 1/5;
sigma = 1/8;
R0=16;

year_table = 1968:2019; 
tau = (R0-.6404)/3.714; %transmission rate is scaled to give R0 of 16 

% alpha = dummy_alpha;
alpha = .1;

year = year_table(1);
filename=['MixingData/ClassMixingData',num2str(year)]
tb1=load(filename,'E1','NGrid', 'tickGrid','ClassProb','DemGrid',...
    'TB','kB','TL','kL','kV','TV','TD','kR','Distrib_Children','StopProb','Vacc_rate');
if(isnan(tb1.Vacc_rate));   Vacc_rate=0; else;  Vacc_rate = tb1.Vacc_rate;    end
tb2=load('ContactMixingData','D_All','D_Ext');
d_int=sum(tb1.ClassProb*(tb2.D_All-tb2.D_Ext));
beta = tau*d_int; % tau is Unit time transmission rate based on Hope-Simpson
Inf_Ext=tau*tb2.D_Ext;

% Set up differential equation to solve for H
dt_H=1;    dQ=10;  run_time=365*1;
options=odeset('RelTol', 1e-3);

%create a struct for storing all values
t =zeros(1,run_time*length(year_table)); 
S_bar=zeros(1,run_time*length(year_table));     E_bar=zeros(run_time*length(year_table),1);
I_bar=zeros(1,run_time*length(year_table));     R_bar=zeros(run_time*length(year_table),1);
yearly_infection=zeros(length(year_table),1);    daily_new_inf=zeros(run_time*length(year_table),1);
H = [];     P_R_yearly=zeros(1,length(year_table));  

[TR,~,~,VectN,Vect,nVect,nTicker,nVectN,E1,E2] = Get_Eq_Demography_seir(popH,tb1.TB,tb1.kB,tb1.TL,tb1.kL,tb1.kV,tb1.TV,tb1.TD,tb1.kR,tb1.Distrib_Children,tb1.StopProb,Vacc_rate,tb1.NGrid,tb1.tickGrid,tb1.E1);

filename = ['alpha/temp_' num2str(year) '_alpha_' num2str(alpha) '.mat']
tb=load(filename,'H_init1','P_R1');
H_init=round(tb.H_init1);
nVectS = nVect(1,:);    nVectE = nVect(2,:);    nVectI = nVect(3,:); nVectR = nVect(4,:);
H_init(nVectI==1 & nVectS==2 & nVectE==0 & nVectR==0 & nTicker==4)=H_init(nVectI==1 & nVectS==2 & nVectE==0 & nVectR==0 & nTicker==4)+1;

[t1,I_bar1,S_bar1,E_bar1,R_bar1,H_init1,P_R1,new_inf] = Get_I_seir_stochastic(run_time,H_init,tb.P_R1,options,dt_H,dQ,VectN,Vect,nVect,nVectN,tb1.NGrid,nTicker,tb1.tickGrid,...
    tb1.StopProb,E1,E2,tb1.DemGrid,alpha,beta,Inf_Ext,g,sigma,Vacc_rate,tb1.kB,tb1.kV,tb1.kL,tb1.kR,tb1.TB,tb1.TV,tb1.TL,TR,seed);

t(1:run_time)=t1;   
S_bar(1:run_time)=S_bar1;    E_bar(1:run_time)=E_bar1;   
I_bar(1:run_time)=I_bar1;    R_bar(1:run_time)=R_bar1;
daily_new_inf(1:run_time)=new_inf; 
cum_inf = cumsum(new_inf);      
yearly_infection(1) = round((cum_inf(end)/sum(round(nVectN.*H_init1')))*1e5);
H{1} = H_init1;     P_R_yearly(1)=P_R1; 

for i=2:length(year_table)
    year = year_table(i);

    filename=['MixingData/ClassMixingData',num2str(year)]
    tb1=load(filename,'E1','NGrid', 'tickGrid','ClassProb','DemGrid',...
        'TB','kB','TL','kL','kV','TV','TD','kR','Distrib_Children','StopProb','Vacc_rate');
    if(isnan(tb1.Vacc_rate));   Vacc_rate=0; else;  Vacc_rate = tb1.Vacc_rate;    end

    maxN = find(tb1.StopProb==1,1);
    H_init=round(H_init1(nVectN<=maxN-1 | (tb1.kB<nTicker & nTicker<=tb1.kB+tb1.kL+tb1.kV & nVectN==maxN)));
    nVectS = nVect(1,:);    nVectE = nVect(2,:);    nVectI = nVect(3,:); nVectR = nVect(4,:);
    H_init(nVectI==1 & nVectS==2 & nVectE==0 & nVectR==0 & nTicker==4)=H_init(nVectI==1 & nVectS==2 & nVectE==0 & nVectR==0 & nTicker==4)+1;
    [TR,~,~,VectN,Vect,nVect,nTicker,nVectN,E1,E2] = Get_Eq_Demography_seir(popH,tb1.TB,tb1.kB,tb1.TL,tb1.kL,tb1.kV,tb1.TV,tb1.TD,tb1.kR,tb1.Distrib_Children,tb1.StopProb,Vacc_rate,tb1.NGrid,tb1.tickGrid,tb1.E1);
    d_int=sum(tb1.ClassProb*(tb2.D_All-tb2.D_Ext));
    beta = tau*d_int; % tau is Unit time transmission rate 

    [t_i,I_bar_i,S_bar_i,E_bar_i,R_bar_i,H_init1,P_R1,new_inf] = Get_I_seir_stochastic(run_time,H_init,P_R1,options,dt_H,dQ,VectN,Vect,nVect,nVectN,tb1.NGrid,nTicker,tb1.tickGrid,...
        tb1.StopProb,E1,E2,tb1.DemGrid,alpha,beta,Inf_Ext,g,sigma,Vacc_rate,tb1.kB,tb1.kV,tb1.kL,tb1.kR,tb1.TB,tb1.TV,tb1.TL,TR,seed);

    t((i-1)*run_time+1:i*run_time)=t_i+t((i-1)*run_time);
    S_bar((i-1)*run_time+1:i*run_time)=S_bar_i;
    E_bar((i-1)*run_time+1:i*run_time)=E_bar_i;
    I_bar((i-1)*run_time+1:i*run_time)=I_bar_i;
    R_bar((i-1)*run_time+1:i*run_time)=R_bar_i;
    daily_new_inf((i-1)*run_time+1:i*run_time)=new_inf;
    cum_inf = cumsum(new_inf);
    yearly_infection(i) = round((cum_inf(end)/sum(round(nVectN.*H_init1')))*1e5);
    H{i} = H_init1;     P_R_yearly(i) = P_R1;
end

filename = ['sim_alpha_' num2str(alpha) '_' datestr(now,'ddmmyy_HHMMSS') '.mat']
save(filename,'t','I_bar','S_bar','R_bar','E_bar','daily_new_inf','yearly_infection','popH','seed','R0','tau','g','sigma','alpha','H','P_R_yearly')
end



