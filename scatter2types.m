% Code for producing scatter plots for graph G

% n=10000;


T=500;
c=5;
p=c/n;
gamma=1;
% Beta=0.1;
eps=0.01;
C=20;
its=5000;
res=1;
fin=2;
% strt=0.09;
strt=0;
% fin=0.85;
fin=fin-strt;


rec1=1;
rho=1;

% G = sparse(rand(n,n)) < p;
% G = sparse(triu(G,1));
% G = G + G';

G2=sparse(2*n);


% figure();
% plot(graph(G));
% figure();
% plot(graph(G2));

[ei,ej]=find(G);
nedge=length(ei);
t11=zeros(nedge,1);
t12=zeros(nedge,1);
t21=zeros(nedge,1);
t22=zeros(nedge,1);



snum=1;

sXa=1:snum;
sXb=1:snum;



grph1=graph(G);

M=[G,speye(n)-spdiags(degree(grph1),0,n,n);speye(n),sparse(n,n)];
[v,eg]=eigs(M,1); v=abs(v);
nbcent=v(1:n);

%        avcent=mean(nbcent);
%        gstrt=find(nbcent>avcent);
scent=sort(nbcent);
gstrt=find(nbcent>scent(floor(699*n/700)));
lgs=length(gstrt);
%%
eigns=zeros(res,1);

osize=zeros(its,1);

thresh=n/700;

nodes=[1:n];
nnodes=length(nodes);
times1=zeros(nnodes,1);
times2=zeros(nnodes,1);
% times3=zeros(nnodes,1);

Vfa=zeros(nnodes,res);
Vfb=zeros(nnodes,res);
% Vfc=zeros(nnodes,res);
%%
for b=1:res
    
    
    G11=sparse(n);
    G21=sparse(n);
    G12=sparse(n);
    G22=sparse(n);
    
    
    Va=zeros(nnodes,1);
    Vb=zeros(nnodes,1);
    
    beta=0.052;
    
    
    Paa=9/10*beta;
    Pab=1/10*beta;
    Pba=3*beta;
    Pbb=1*beta;
    
    Tmat=[Paa/(Paa+Pab+rho),Pba/(Pba+Pbb+rho);Pab/(Paa+Pab+rho),Pbb/(Pba+Pbb+rho)];
    [E,e]=eigs(Tmat,1);
    E=abs(E);
    eigns(b)=e;
    
    
    for i=1:its
        
        
        i/its
        
        snode=gstrt(randi(lgs));
        stype=    randsample(2,1,true,E)-1;
        snode=snode+n*stype;
        
        t11=exprnd(1/Paa,nedge,1);
        t12=exprnd(1/Pab,nedge,1);
        t21=exprnd(1/Pba,nedge,1);
        t22=exprnd(1/Pbb,nedge,1);
        r=exprnd(1,nedge,1);
        %
        gt11=t11;
        gt12=t12;
        gt21=t21;
        gt22=t22;
        %      ----- exc ------
        gt11(t11>min(r,t12))=0;
        gt12(t12>min(r,t11))=0;
        gt21(t21>min(r,t22))=0;
        gt22(t22>min(r,t21))=0;
        %    ------ nexc-----
        %                 gt11(t11>r)=0;
        %             gt12(t12>r)=0;
        %             gt21(t21>r)=0;
        %             gt22(t22>r)=0;
        %
        G11=sparse(ei,ej,gt11,n,n);
        G12=sparse(ei,ej,gt12,n,n);
        G21=sparse(ei,ej,gt21,n,n);
        G22=sparse(ei,ej,gt22,n,n);
        
        
        
        
        
        G2=[G11,G12;G21,G22];
        
        grph=digraph(G2);
        
        [~,times]= shortestpathtree(grph,snode,1:2*n);
        times1=times(1:n);
        times2=times(n+1:2*n);
        
        % %   toc;
        %  ------  exc-------
        gt1=find(times1<times2);
        gt2=find(times2<times1);
        % --------   nexc --------
        %               gt1=find(times1~=Inf);
        %             gt2=find(times2~=Inf);
        
        
        osize(i)=length(gt1)+length(gt2);
        if osize(i)>thresh
            Va(gt1)=Va(gt1)+1;
            Vb(gt2)=Vb(gt2)+1;
        end
        
        
    end
    
    
    Vfa(:,b)=Va./its;
    Vfb(:,b)=Vb./its;
    
    
end

%%




figure();hold on;

plot(nbcent,Vfb,'r.','LineWidth',1);
plot(nbcent,Vfa,'b.','LineWidth',1);
xlabel('nb cent');
ylabel('Infection probability');