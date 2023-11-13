function [R] = convert_tau_to_R0

%params fixed
popH = 5e4; %dummy_pop
g = 1/5;
sigma = 1/8;

year = 2000;
tau = 4;

filename=['MixingData/ClassMixingData',num2str(year)];
load(filename,'ClassProb', 'E1','NGrid', 'tickGrid','DemGrid',...
    'TB','kB','TL','kL','kV','TV','TD','kR','Distrib_Children','StopProb');
load('ContactMixingData','D_All','D_Ext');

Vacc_rate = 0;
maxN = find(StopProb==1,1);

d_int=sum(ClassProb*(D_All-D_Ext));
beta = tau*d_int; % tau is Unit time transmission rate based on Hope-Simpson
Inf_Ext=tau*D_Ext;
m_run = maxN*(1/g+1/sigma);
r = Get_r0_seir(Inf_Ext,beta,g,sigma,popH,kB,kL,kV,kR,TB,TL,TV,TD,Distrib_Children,...
    StopProb,Vacc_rate,DemGrid,NGrid,tickGrid,E1,m_run,0);
R =  (1+r/g)*(1+r/sigma);

end 
