function [Q_demo]=Get_Qdemo_seir(P_R,VectN,Vect,Vacc_rate,kB,kV,kL,kR,TB,TV,TL,TR,StopProb)

VectS=Vect(1,:);    VectE=Vect(2,:);    VectI=Vect(3,:);    VectR=Vect(4,:);    maxN = max(VectN);
        
%Set up rates
Birth_Ticker_Rates=ones(maxN,kB)*kB/TB;
Prob_at_Birth=sparse(1+kB+0*StopProb, 1:maxN, StopProb, 1+kB, maxN);
Prob_at_Birth=Prob_at_Birth+sparse(1+0*StopProb, 1:maxN, 1-StopProb, 1+kB,maxN);
Prob_at_Birth(1,maxN)=0; Prob_at_Birth(1+kB,maxN)=1;

Vacc_Ticker_Rates=ones(maxN,kV)*kV/TV;
Leave_Ticker_Rates=ones(maxN,kL)*kL./(TL - ([1:maxN]'-3)*ones(1,kL)*TB);
% Expected intervals for subsequent leavings are just birth intervals
Leave_Ticker_Rates2=ones(maxN,kB)*kB./(TB);
Leave_Ticker_Rates2(end,:)=0; % as cannot be in state maxN as one child has left

Reset_Ticker_Rates=ones(maxN,kR)*kR./(TR);
%Set up what to do after a reset
if max(VectR)<2 % Inf_at_Reset is probability of a child having been infected before starting new home
    Inf_at_Reset=0;
else
    Inf_at_Reset=P_R; % Proportion of newly recovered households
end

Tickers=(kB)+(kV)+(kL)+(kB)+(kR); % Total number of possible tick values
L=size(Vect,2); SIZE=Tickers*L; % L is number of infectious states in max size household

for i=Tickers:-1:1 % Puts in full range of ticker classes for each epidemic state
    j=i-1;
    nVectN(j*L+[1:L])=VectN;
    nTicker(j*L+[1:L])=i;
end

Q_demo=sparse(1,1,0,SIZE,SIZE); % Define demographic transition matrix

% Now do Births & Ticker
for i=1:(kB-1) % Again, this is for each k class to up-date the ticker
    j=i-1;
    %first remove the rate of advancement from the diagonal.
    From=[1:L]+j*L; To=[1:L]+j*L; Rate=-Birth_Ticker_Rates(VectN(1:L),i)';
    pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);

    From = [1:L];
    for k=1:length(From)
        Rate = Birth_Ticker_Rates(VectN(From(k)),i);

        if VectN(From(k))==2 % no kids
            Q_demo(From(k)+j*L,From(k)+j*L+L)=Q_demo(From(k)+j*L,From(k)+j*L+L)+Rate;
        elseif Vect(1,From(k)) == 0 % no susceptible in the household
            Q_demo(From(k)+j*L,From(k)+j*L+L)=Q_demo(From(k)+j*L,From(k)+j*L+L)+Rate;
        else % move the kid (one from S to R)
            tmpV = Vect(:,From(k)); tmpV = tmpV+[-1;0;0;1];
            m=find(sum(abs(Vect-tmpV*ones(1,L)),1)<0.01);
            To = m;
            Q_demo(From(k)+j*L,To+j*L+L)=Q_demo(From(k)+j*L,To+j*L+L)+Rate*Vacc_rate;
            Q_demo(From(k)+j*L,From(k)+j*L+L)=Q_demo(From(k)+j*L,From(k)+j*L+L)+Rate*(1-Vacc_rate);
        end
    end
end

% NOW DO IT FOR BIRTHS
i=kB; j=i-1; 
% remove the diagonals.
From=[1:L]+j*L; To=[1:L]+j*L; Rate=-Birth_Ticker_Rates(VectN(1:L),i)';
pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
From=[1:L]; 
for k=1:length(From)
    tmpVS=Vect(:,From(k)); tmpVS(1)=tmpVS(1)+1; %Everyone susceptible at birth
    Rate=Birth_Ticker_Rates(VectN(From(k)),i);
    Stop_Prob=Prob_at_Birth(end,VectN(From(k)));    % This is simply referencing StopProb vector
    
    if VectN(From(k))==2 % stopping before any kids
        To=From(k); t=(kB)+(kV)+(kL)+(kB);
        Q_demo(From(k)+j*L,To+t*L)=Q_demo(From(k)+j*L,To+t*L)+Rate*Stop_Prob;
    else % stopping after this kid  -- since N is greater than 2, they have already 
        % passed through one kB, hence t=kB for stopping criteria
        To=From(k); t=(kB);
        Q_demo(From(k)+j*L,To+t*L)=Q_demo(From(k)+j*L,To+t*L)+Rate*Stop_Prob;
    end
    
    % Adding a child
    if Stop_Prob<1
        m=find(sum(abs(Vect-tmpVS*ones(1,L)),1)<0.01); % adding a susceptible
        To=m; t=0; %reset ticker
        if VectN(To)==maxN
            t=kB;
        end
        Q_demo(From(k)+j*L,To+t*L)=Q_demo(From(k)+j*L,To+t*L)+Rate*(1-Stop_Prob);
        
    end
end
        
        
% LAST CHILD VACCINATING
% First ticker increase after all the children are born     % Important step, made an error here previously

F=find(VectN>2);
i=kV;   j=i+kB-1;    %j=kB
From=F+j*L; To=F+j*L; Rate=-Vacc_Ticker_Rates(VectN(F),i)';
pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);

for k = 1:length(F)
    Rate = Vacc_Ticker_Rates(VectN(F(k)),i)'; 
    if VectS(F(k)) > 0 
        tmpV = Vect(:,F(k)); tmpV = tmpV+[-1;0;0;1];
        m=find(sum(abs(Vect-tmpV*ones(1,L)),1)<0.01);
        To = m;
        Q_demo(F(k)+j*L,To+j*L+L)=Q_demo(F(k)+j*L,To+j*L+L)+ Rate*Vacc_rate;
        Q_demo(F(k)+j*L,F(k)+j*L+L)=Q_demo(F(k)+j*L,F(k)+j*L+L)+ Rate*(1-Vacc_rate);
    else
        Q_demo(F(k)+j*L,F(k)+j*L+L)=Q_demo(F(k)+j*L,F(k)+j*L+L)+ Rate;
    end 
end 

% Now do First Leaving & Ticker
% FIRST CHILD LEAVING
F=find(VectN>2);
for i=1:(kL-1)
    j=i+kB+kV-1;
    From=F+j*L; To=F+j*L; Rate=-Leave_Ticker_Rates(VectN(F),i)';
    pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
    
    From=F+j*L; To=F+j*L+L; Rate=Leave_Ticker_Rates(VectN(F),i)';
    pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
end

i=kL; j=i+kB+kV-1;
%Now do what happens at leaving
From=find(VectN>2); clear To;
for k=1:length(From)
      Rate=Leave_Ticker_Rates(VectN(From(k)),i);
      tmpV=Vect(:,From(k)); tmpN = sum(tmpV);
      
      RateE = tmpV(2)/tmpN;
      RateI = tmpV(3)/tmpN;
      RateR = tmpV(4)/tmpN;
%       if tmpV(4)==0
%           RateR = 0;
%           RateE = tmpV(2)/tmpN;
%           RateI = tmpV(3)/tmpN;
%       elseif tmpV(4)==1
%           RateR = (1-P_R)/(tmpN-P_R);
%           RateE = tmpV(2)/(tmpN-P_R);
%           RateI = tmpV(3)/(tmpN-P_R);
%       else
%           RateR = (tmpV(4)-2*P_R)/(tmpN-2*P_R);
%           RateE = tmpV(2)/(tmpN-2*P_R);
%           RateI = tmpV(3)/(tmpN-2*P_R);
%       end
      RateS = (1-RateE-RateI-RateR);

      if VectS(From(k))>0
          tmpVS = tmpV; tmpVS(1) = tmpVS(1)-1;  % losing a susceptible
          m=find(sum(abs(Vect-tmpVS*ones(1,L)),1)<0.01);
          ToS = m(1);
          if VectN(From(k))==3
              t=(kB)+(kV)+(kL)+(kB);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToS+t*L,RateS*Rate,SIZE,SIZE); 
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateS*Rate,SIZE,SIZE);
          else
              t=(kB)+(kV)+(kL);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToS+t*L,RateS*Rate,SIZE,SIZE); 
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateS*Rate,SIZE,SIZE);
          end
      end
      if VectE(From(k))>0
          tmpVE = tmpV; tmpVE(2) = tmpVE(2)-1;      % losing an E class indv
          m=find(sum(abs(Vect-tmpVE*ones(1,L)),1)<0.01);
          ToE = m(1);
          if VectN(From(k))==3
              t=(kB)+(kV)+(kL)+(kB);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToE+t*L,RateE*Rate,SIZE,SIZE); 
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateE*Rate,SIZE,SIZE);
          else
              t=(kB)+(kV)+(kL);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToE+t*L,RateE*Rate,SIZE,SIZE); 
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateE*Rate,SIZE,SIZE);
          end
      end
      if VectI(From(k))>0
          tmpVI = tmpV; tmpVI(3) = tmpVI(3)-1;      % losing an infected
          m=find(sum(abs(Vect-tmpVI*ones(1,L)),1)<0.01);
          ToI = m(1);
          if VectN(From(k))==3
              t=(kB)+(kV)+(kL)+(kB);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToI+t*L,RateI*Rate,SIZE,SIZE); 
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateI*Rate,SIZE,SIZE);
          else
              t=(kB)+(kV)+(kL);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToI+t*L,RateI*Rate,SIZE,SIZE); 
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateI*Rate,SIZE,SIZE);
          end
      end
      if VectR(From(k))>0
          tmpVR = tmpV; tmpVR(4) = tmpVR(4)-1;      % losing a recovered
          m=find(sum(abs(Vect-tmpVR*ones(1,L)),1)<0.01);
          ToR = m(1);
          if VectN(From(k))==3
              t=(kB)+(kV)+(kL)+(kB);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToR+t*L,RateR*Rate,SIZE,SIZE);
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateR*Rate,SIZE,SIZE);
          else
              t=(kB)+(kV)+(kL);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToR+t*L,RateR*Rate,SIZE,SIZE);
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateR*Rate,SIZE,SIZE);
          end
      end
      
end
   
 % NOW OTHER CHILDREN LEAVING     
               
F=find(VectN>2);
for i=1:(kB-1)
    j=i+kB+kV+kL-1;
    From=F+j*L; To=F+j*L; Rate=-Leave_Ticker_Rates2(VectN(F),i)';
    pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
    
    From=F+j*L; To=F+j*L+L; Rate=Leave_Ticker_Rates2(VectN(F),i)';
    pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
end

i=kB; j=i+kB+kV+kL-1;
%Now do what happens at leaving
From=find(VectN>2); clear To;
for k=1:length(From)
      Rate=Leave_Ticker_Rates2(VectN(From(k)),i);
      tmpV=Vect(:,From(k)); tmpN = sum(tmpV);

      RateE = tmpV(2)/tmpN;
      RateI = tmpV(3)/tmpN;
      RateR = tmpV(4)/tmpN;
%       if tmpV(4)==0
%           RateR = 0;
%           RateE = tmpV(2)/tmpN;
%           RateI = tmpV(3)/tmpN;
%       elseif tmpV(4)==1
%           RateR = (1-P_R)/(tmpN-P_R);
%           RateE = tmpV(2)/(tmpN-P_R);
%           RateI = tmpV(3)/(tmpN-P_R);
%       else
%           RateR = (tmpV(4)-2*P_R)/(tmpN-2*P_R);
%           RateE = tmpV(2)/(tmpN-2*P_R);
%           RateI = tmpV(3)/(tmpN-2*P_R);
%       end
      RateS = (1-RateE-RateI-RateR);
      
      if VectS(From(k))>0
          tmpVS = tmpV; tmpVS(1) = tmpVS(1)-1;         % losing a susceptible
          m=find(sum(abs(Vect-tmpVS*ones(1,L)),1)<0.01);
          ToS = m(1);
          if VectN(From(k))==3
              t=(kB)+(kV)+(kL)+(kB);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToS+t*L,RateS*Rate,SIZE,SIZE); 
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateS*Rate,SIZE,SIZE);
          else
              t=(kB)+(kV)+(kL);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToS+t*L,RateS*Rate,SIZE,SIZE); 
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateS*Rate,SIZE,SIZE);
          end
      end
      if VectE(From(k))>0
          tmpVE = tmpV; tmpVE(2) = tmpVE(2)-1;      % losing an E class indv
          m=find(sum(abs(Vect-tmpVE*ones(1,L)),1)<0.01);
          ToE = m(1);
          if VectN(From(k))==3
              t=(kB)+(kV)+(kL)+(kB);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToE+t*L,RateE*Rate,SIZE,SIZE); 
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateE*Rate,SIZE,SIZE);
          else
              t=(kB)+(kV)+(kL);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToE+t*L,RateE*Rate,SIZE,SIZE); 
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateE*Rate,SIZE,SIZE);
          end
      end
      if VectI(From(k))>0
          tmpVI = tmpV; tmpVI(3) = tmpVI(3)-1;       % losing an infected
          m=find(sum(abs(Vect-tmpVI*ones(1,L)),1)<0.01);
          ToI = m(1);
          if VectN(From(k))==3
              t=(kB)+(kV)+(kL)+(kB);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToI+t*L,RateI*Rate,SIZE,SIZE); 
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateI*Rate,SIZE,SIZE);
          else
              t=(kB)+(kV)+(kL);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToI+t*L,RateI*Rate,SIZE,SIZE); 
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateI*Rate,SIZE,SIZE);
          end
      end
      if VectR(From(k))>0
          tmpVR = tmpV; tmpVR(4) = tmpVR(4)-1;       % losing a recovered
          m=find(sum(abs(Vect-tmpVR*ones(1,L)),1)<0.01);
          ToR = m(1);
          if VectN(From(k))==3
              t=(kB)+(kV)+(kL)+(kB);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToR+t*L,RateR*Rate,SIZE,SIZE);
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateR*Rate,SIZE,SIZE);
          else
              t=(kB)+(kV)+(kL);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToR+t*L,RateR*Rate,SIZE,SIZE);
              Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateR*Rate,SIZE,SIZE);
          end
      end
end

       
% NOW DO RESET AND TICKER
F=find(VectN==2); % must be just 2 adults.
for i=1:(kR-1)
    j=i+kB+kV+kL+kB-1;
    From=F+j*L; To=F+j*L; Rate=-Reset_Ticker_Rates(VectN(F),i)';
    pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
    
    From=F+j*L; To=F+j*L+L; Rate=Reset_Ticker_Rates(VectN(F),i)';
    pos=find(From>0 & To>0 & Rate~=0); 
    Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
end

%Now do what happens at reset
i=kR; j=i+kB+kV+kL+kB-1;
From=find(VectN==2); clear To;
ToProb(3)=(1-Inf_at_Reset).^2; To(3)=0; % both sus
m=find(VectS==2 & VectE==0 & VectI==0 & VectR==0); % Passed over by outbreak
if ~isempty(m)
    To(3)=m;
end

ToProb(2)=2*(Inf_at_Reset)*(1-Inf_at_Reset);
m=find(VectS==1 & VectE==0 & VectI==0 & VectR==1); % 1 recovered
if ~isempty(m)
    To(2)=m;
end

ToProb(1)=(Inf_at_Reset).^2;
m=find(VectS==0 & VectE==0 & VectI==0 & VectR==2); % 2 recovers
if ~isempty(m)
    To(1)=m;
end

Rate=Reset_Ticker_Rates(VectN(From),i);
for k=1:length(From)
    Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-Rate(k),SIZE,SIZE);
    for K=1:3
       if Rate(k)*ToProb(K)>0 
           Q_demo=Q_demo+sparse(From(k)+j*L,To(K),Rate(k)*ToProb(K),SIZE,SIZE);% N=N+sparse(From(k),From(k),-Rate(k)*ToProb(K),SIZE,SIZE);
       end
    end
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

Q_demo(r,:)=[]; Q_demo(:,r)=[];
end 


