
%Non exclusive infection with 3 types - simulation

n=1000;
c=5;
p=c/n;


its=1500;
res=50;
fin=0.15;
strt=0;
fin=fin-strt;

rho=1;

% G = sparse(rand(n,n)) < p;
% G = sparse(triu(G,1));
% G = G + G';

G2=sparse(2*n);



[ei,ej]=find(G);
nedge=length(ei);
t11=zeros(nedge,1);
t12=zeros(nedge,1);
t21=zeros(nedge,1);
t22=zeros(nedge,1);

times1=zeros(n,1);
times2=zeros(n,1);
times3=zeros(n,1);


Vfa=zeros(n,res,its);
Vfb=zeros(n,res,its);
Vfc=zeros(n,res,its);


%%
eigns=zeros(res,1);


parfor b=1:res
    %%
    
    G11=sparse(n);
G21=sparse(n);
G12=sparse(n);
G22=sparse(n);
    Va=zeros(n,its);
    Vb=zeros(n,its);
    Vc=zeros(n,its);

    
    beta=exp(b/10)*fin/exp(res/10)+strt



p11=8/10*beta;
     p12=1/10*beta;
     p13=1/10*beta;
     p21=2*beta;
     p22=6*beta;
     p23=2*beta;
     p31=10*beta;
     p32=10*beta;
     p33=10*beta;
     
Tmat=[p11/(p11+rho),p21/(p21+rho),p31/(p31+rho);p12/(p12+rho),p22/(p22+rho),p32/(p32+rho);p13/(p13+rho),p23/(p23+rho),p33/(p33+rho)];


      [E,e]=eigs(Tmat,1);
     EE= sum(E);
      E=E./EE; 
       

      eigns(b)=e;
     

for i=1:its
    %%
   
    snode=randi(n);
stype=    randsample(3,1,true,E)-1;
snode=snode+n*stype;

    t11=exprnd(1/p11,nedge,1);
     t12=exprnd(1/p12,nedge,1);
     t13=exprnd(1/p13,nedge,1);
      t21=exprnd(1/p21,nedge,1);
     t22=exprnd(1/p22,nedge,1);
     t23=exprnd(1/p23,nedge,1);
       t31=exprnd(1/p31,nedge,1);
     t32=exprnd(1/p32,nedge,1);
     t33=exprnd(1/p33,nedge,1);
     %%%%r=exprnd(1,nedge,1);
     %rnode=exprnd(1,n,1);
     %r=zeros(nedge,1);
     %for node=1:n
     %    r(find(ei==node))=rnode(node);
     %end
     r1node=exprnd(1,n,1);
     r1=zeros(nedge,1);
     for node=1:n
         r1(find(ei==node))=r1node(node);
     end   
     r2node=exprnd(1,n,1);
     r2=zeros(nedge,1);
     for node=1:n
         r2(find(ei==node))=r2node(node);
     end  
     r3node=exprnd(1,n,1);
     r3=zeros(nedge,1);
     for node=1:n
         r3(find(ei==node))=r3node(node);
     end  
     %     

     gt11=t11;
     gt12=t12;
     gt13=t13;
     gt21=t21;
     gt22=t22;
     gt23=t23;
     gt31=t31;
     gt32=t32;
     gt33=t33;
     
%      ----- exc ------
%     gt11(t11>min(r,t12))=0;
%     gt12(t12>min(r,t11))=0;
%     gt21(t21>min(r,t22))=0;
%     gt22(t22>min(r,t21))=0;
%    ------ nexc-----
     gt11(t11>r1)=0;
    gt12(t12>r1)=0;
     gt13(t13>r1)=0;
    gt21(t21>r2)=0;
    gt22(t22>r2)=0;
     gt23(t23>r2)=0;
    gt31(t31>r3)=0;
    gt32(t32>r3)=0;
     gt33(t33>r3)=0;
%     
    G11=sparse(ei,ej,gt11,n,n);
    G12=sparse(ei,ej,gt12,n,n);
    G13=sparse(ei,ej,gt13,n,n);
    G21=sparse(ei,ej,gt21,n,n);
    G22=sparse(ei,ej,gt22,n,n);
    G23=sparse(ei,ej,gt23,n,n);
    G31=sparse(ei,ej,gt31,n,n);
    G32=sparse(ei,ej,gt32,n,n);
    G33=sparse(ei,ej,gt33,n,n);

    
    G2=[G11,G12,G13;G21,G22,G23;G31,G32,G33];
    
    grph=digraph(G2);
    

   [~,times]= shortestpathtree(grph,snode,1:3*n);
  times1=times(1:n);
    times2=times(n+1:2*n);
    times3=times(2*n+1:3*n);

%  ------  exc-------
%     gt1=find(times1<times2);
%     gt2=find(times2<times1);
% --------   nexc --------
      gt1=find(times1<1000);
    gt2=find(times2<1000);
gt3=find(times3<1000);


   Va(gt1,i)= 1;
            Vb(gt2,i)= 1;
             Vc(gt3,i)= 1;
          
end

 Vfa(:,b,:)=Va;
    Vfb(:,b,:)=Vb;
     Vfc(:,b,:)=Vc;
end

%%

thresh=50;
osize=sum(Vfa,1)+sum(Vfb,1)+sum(Vfc,1);
shift=0.0;

Vfa2=zeros(n,res);
Vfb2=zeros(n,res);
Vfc2=zeros(n,res);
for b=1:res
    gits=nnz(osize(:,b,:)>thresh);
    if gits>its/10
Vfa2(:,b)=sum(Vfa(:,b,osize(:,b,:)>thresh),3)./gits;
Vfb2(:,b)=sum(Vfb(:,b,osize(:,b,:)>thresh),3)./gits;
Vfc2(:,b)=sum(Vfc(:,b,osize(:,b,:)>thresh),3)./gits;
    end
end

% 
 degs=sum(G);

Vfa3=Vfa2(degnodes,:);
Vfb3=Vfb2(degnodes,:);
Vfc3=Vfc2(degnodes,:);
 %%
for i=1:2
    frst=find(Vfa3(i,:)>0.001,1);
Vfa3(i,frst:res)=smooth(Vfa3(i,frst:res),1);
Vfb3(i,frst:res)=smooth(Vfb3(i,frst:res),1);
Vfc3(i,frst:res)=smooth(Vfc3(i,frst:res),1);
end
figure();hold on;
plot(eigns(frst-1:end),Vfb3(:,frst-1:end),'r','LineWidth',1);
 plot(eigns(frst-1:end),Vfa3(:,frst-1:end),'b','LineWidth',1);
  plot(eigns(frst-1:end),Vfc3(:,frst-1:end),'g','LineWidth',1);
   xlabel('$\lambda_T$');
ylabel('Infection probability');
axis([0,1,0,1]);

% figure();hold on;
%  plot(eigns,Vfa3+Vfb3+Vfc3,'color',[0.8,0.5,0.8],'LineWidth',1);
%    xlabel('\lambda_T');
% ylabel('Infection probability');