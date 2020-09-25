function [gt1,gt2] = simrun(bet1,bet2,mu1,mu2,rho,ei,ej,nedge,n)

%simulation run for phase diagram

 Paa=(1-mu1)*bet1;
    Pab=mu1*bet1;
    Pbb=(1-mu2)*bet2;
    Pba=mu2*bet2;

Tmat=[Paa/(Paa+Pab+rho),Pba/(Pba+Pbb+rho);Pab/(Paa+Pab+rho),Pbb/(Pba+Pbb+rho)];
    [E,~]=eigs(Tmat,1);
     EE= sum(E);
      E=E./EE; 

      
snode=randi(n);
stype=    randsample(2,1,true,E)-1;
snode=snode+n*stype;

    t11=exprnd(1/Paa,nedge,1);
     t12=exprnd(1/Pab,nedge,1);
      t21=exprnd(1/Pba,nedge,1);
     t22=exprnd(1/Pbb,nedge,1);
     r=exprnd(1/rho,nedge,1);
    
  
     
    t11(t11>min(r,t12))=0;
    t12(t12>min(r,t11))=0;
    t21(t21>min(r,t22))=0;
    t22(t22>min(r,t21))=0;
    
    G11=sparse(ei,ej,t11,n,n);
    G12=sparse(ei,ej,t12,n,n);
    G21=sparse(ei,ej,t21,n,n);
    G22=sparse(ei,ej,t22,n,n);
    

    
    G2=[G11,G12;G21,G22];
    
    grph=digraph(G2);
    
  

   [~,times]= shortestpathtree(grph,snode,1:2*n);
  times1=times(1:n);
    times2=times(n+1:2*n);

   
    gt1=nnz(times1<times2);
    gt2=nnz(times2<times1);
    


end

