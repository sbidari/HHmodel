function [ClassProb,E1,NGrid,tickGrid,DemGrid,TB,kB,TL,kL,kV,TV,TD,kR,Distrib_Children,StopProb,Vacc_rate] = ClassMixing

ModelBds=[0 1 14 50 99];

%Use world data
Data = readtable('dataInput/world.csv','ReadRowNames',true);

year = 2015

TB=2*365;
kB=2;
TV=1*365;
kV=1;
TL=Data{num2str(year),'MarriageAge'}*365-TV;
kL=1;
TD=Data{num2str(year),'LifeExpectancyAtBirth'}*365;
kR=2;
DC=Data{num2str(year),'ChildProb'};
Vacc_rate = Data{num2str(year),'VC'}/100;

Distrib_Children=str2num(cell2mat(DC));

if (year>=1993) 
    Constrained_max=5; % This is the most kids we can have under our model assumptions with these parameter choices
    if length(Distrib_Children)>Constrained_max+1
        Distrib_Children(Constrained_max+1)=sum(Distrib_Children(Constrained_max+1:end)); % +1 since first one is no kids
        Distrib_Children=Distrib_Children(1:Constrained_max+1)
    end
elseif (year<1993) &&  (year >1977)
    Constrained_max=6; % This is the most kids we can have under our model assumptions with these parameter choices
    if length(Distrib_Children)>Constrained_max+1
        Distrib_Children(Constrained_max+1)=sum(Distrib_Children(Constrained_max+1:end)); % +1 since first one is no kids
        Distrib_Children=Distrib_Children(1:Constrained_max+1)
    end
else 
    Constrained_max=7;
    if length(Distrib_Children)>Constrained_max+1
        Distrib_Children(Constrained_max+1)=sum(Distrib_Children(Constrained_max+1:end)); % +1 since first one is no kids
        Distrib_Children=Distrib_Children(1:Constrained_max+1)
    end
end


Distrib_Children=Distrib_Children/sum(Distrib_Children);
Cond_DC=Distrib_Children(2:end)/sum(Distrib_Children(2:end)); % Distribution of number of children conditioned on having kids
tmp=[0 cumsum(Distrib_Children)];
StopProb=[0 Distrib_Children./(1-tmp(1:(end-1)))]; % Prob of stopping after each number of kids
StopProb(end)=1;

disp('Now calculating translation matrix')
points=10000; % This is number of points for Monte Carlo integration
[E1,NGrid,tickGrid,DemGrid] = State_to_Class_Mat_With_Elders(kB,kV,kL,kR,TB,TV,TL,365*ModelBds,StopProb,points);
disp('Translation matrix calculated')

disp('Now calculating age class proportions')
ClassProb=zeros(1,length(ModelBds)-1); % This loop will do age class proportions for children
for kk=1:length(ModelBds)-1
    tic
    L=365*ModelBds(kk); U=365*ModelBds(kk+1);
    fun1 = @(A,L,kB,kL,kV,TB,TL,TV,Cond_DC) LeaveAgeProb(A,kB,kL,kV,TB,TL,TV,Cond_DC).*((A-L)./A); % Integrand for leaving within age class
    fun2 = @(A,L,U,kB,kL,kV,TB,TL,TV,Cond_DC) LeaveAgeProb(A,kB,kL,kV,TB,TL,TV,Cond_DC).*((U-L)./A); % Integrand for leaving after age class
    int1=integral(@(A)fun1(A,L,kB,kL,kV,TB,TL,TV,Cond_DC),L,U);
    int2=integral(@(A)fun2(A,L,U,kB,kL,kV,TB,TL,TV,Cond_DC),U,150*365); % Assume people have definitely left home by LE
    ClassProb(kk)=int1+int2;
    toc
end
% Now assume adults spend TR-TL alive after leaving home
ClassProb=ClassProb*(TL+TV)/TD; % Scale by proportion of lifetime spent as child
ClassProb(end-1)=ClassProb(end-1)+(365*ModelBds(end-1)-(TL+TV))/TD;
ClassProb(end)= 1-sum(ClassProb);
disp('Age class proportions calculated')

filename=['ClassMixingData',num2str(year)]
save(filename,'ClassProb','E1','NGrid','tickGrid','DemGrid','TB','kB','TL','kL','kV','TV','TD','kR','Distrib_Children','StopProb','Vacc_rate');
end


function [P_Age] = LeaveAgeProb(A,kB,kL,kV,TB,TL,TV,Cond_DC)
% Calculates probability of kid leaving home at age A by going over all
% possible numbers of children and positions
% Recall that Cond_DC is conditioned on having at least one child!
Nmax=length(Cond_DC);
P_Age=0;
for N=1:Nmax
    summand=0;
    for n=1:N
        lam=[(kB/TB)*ones(1,(N-n)*kB) (kV/TV)*ones(1,kV) (kL/(TL-(N-1)*TB))*ones(1,kL) (kB/TB)*ones(1,(n-1))*kB];
        summand=summand+(1/N)*HypoExpPdf(A,lam);
    end
    P_Age=P_Age+Cond_DC(N)*summand;
end
end

function [P_Hyp] = HypoExpPdf(x,lambda)
% Hypoexponential distribution - see wikipedia etc
range = length(x);
n = length(lambda);

alpha = zeros(1,n);
alpha(1) = 1;

index = 1:n;
Theta = sparse(index,index,-lambda,n,n);
Theta = Theta + sparse(index(1:end-1),index(1:end-1)+1,lambda(1:end-1),n,n);

for l=range:-1:1
    xTheta = x(l)*Theta;
    xTheta = expm(xTheta);
    P_Hyp(l) = -alpha*xTheta*Theta*ones(n,1);
end
P_Hyp = full(P_Hyp);
end