
%Non exclusive infection with 3 types - prediction

n=1000;
p=5/n;
C=100;
nsteps=50;
strt=0;
fin=0.15;

numt=10;


%G = rand(n,n) < p; 
%G = triu(G,1); 
%G = G + G'; 


hsize=sum(sum(G));
hash=zeros(hsize);


[edgei,edgej]=find(G==1);
edges=[edgei,edgej];

for e=1:length(edges)
    for f=1:length(edges)
        if edges(e,2)==edges(f,1) && edges(e,1)~=edges(f,2)
            hash(e,f)=1;
        end
    end
end


 rho=1;
 
 V1full=zeros(n,nsteps);
 V2full=zeros(n,nsteps);
 V3full=zeros(n,nsteps);
t1spread=zeros(1,nsteps);

 
  Neighs=cell(hsize,1);
     
     for i=1:hsize
         neighs=find(hash(:,i)==1);
         Neighs{i}=neighs;
     end
     
   
     Neighs2=cell(n);
      
     for i=1:n
         ind=find(edgej==i);
         Neighs2{i}=ind;
     end
     

     %%
     
     
      H1=zeros(hsize);
      H2=zeros(hsize);
        H3=zeros(hsize);
     
 for b=nsteps:-1:1
     %%
     b

beta=exp(b/10)*fin/exp(nsteps/10)+strt
     
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


     

     %%
     for k=2:numt+1
        
         for j=1:hsize
             H1(j)=1-Tmat(1,1)*(1-prod(H1(Neighs{j})))-Tmat(1,2)*(1-prod(H2(Neighs{j})))-Tmat(1,3)*(1-prod(H3(Neighs{j})));
             H2(j)=1-Tmat(2,1)*(1-prod(H1(Neighs{j})))-Tmat(2,2)*(1-prod(H2(Neighs{j})))-Tmat(2,3)*(1-prod(H3(Neighs{j})));
             H3(j)=1-Tmat(3,1)*(1-prod(H1(Neighs{j})))-Tmat(3,2)*(1-prod(H2(Neighs{j})))-Tmat(3,3)*(1-prod(H3(Neighs{j})));
        
         end

     end

        V1=zeros(1,n);
     V2=zeros(1,n);
     V3=zeros(1,n);
     

for i=1:n
         V1(i)=(1-prod(H1(Neighs2{i})));
         V2(i)=(1-prod(H2(Neighs2{i})));
          V3(i)=(1-prod(H3(Neighs2{i})));
end



     %%

     %

     
     V1full(:,b)=V1; 
     V2full(:,b)=V2;
     V3full(:,b)=V3;

     t1spread(b)=e;
 end
%%  

 degs=sum(G);

degnodes=zeros(2,1);
for i=1:2
    degnodes(i)=find(degs==i*6-3,1);
end
frst=find(V2full(1,:)>0.001,1);
figure();hold on;
plot(t1spread(frst-1:end),V2full(degnodes,frst-1:end),'r','LineWidth',1);
 plot(t1spread(frst-1:end),V1full(degnodes,frst-1:end),'b','LineWidth',1);
  plot(t1spread(frst-1:end),V3full(degnodes,frst-1:end),'g','LineWidth',1);
    xlabel('$\lambda_T$');
ylabel('Infection probability');
axis([0,1,0,1]);