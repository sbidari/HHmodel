function [DiseaseFree]=Get_initial_inf_seir(nVect,nTicker,DiseaseFree)
% Adds some infection to initialize the system
    nVectS = nVect(1,:);    nVectE = nVect(2,:);    nVectI = nVect(3,:); nVectR = nVect(4,:); nVectN = sum(nVect);
    DiseaseFree(nVectI==1 & nVectS==2 & nVectE==0 & nVectR==0 & nTicker==4)=1;

end
