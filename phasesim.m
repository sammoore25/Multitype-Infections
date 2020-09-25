%phase diagram - simulation


n=1000;


T=500;
c=4;
p=c/n;
rho=1;
mu1=1/10;
mu2=1/3;

its=500;



% G = sparse(rand(n,n)) < p;
% G = sparse(triu(G,1));
% G = G + G';
% [ei,ej]=find(G);
% nedge=length(ei);



fin1=0.6;res1=fin1/20;strt1=0;
fin2=2.2;res2=fin2/20;strt2=0;

M1=zeros((fin1-strt1)/res1,(fin2-strt2)/res2,its);
M2=zeros((fin1-strt1)/res1,(fin2-strt2)/res2,its);
for b1=1:(fin1-strt1)/res1
    parfor b2=1:(fin2-strt2)/res2
        bet1=b1*res1+strt1;
        bet2=b2*res2+strt2;
        
        gt11=0;gt22=0;

        for i=1:its
            [gt11,gt22]=simrun(bet1,bet2,mu1,mu2,rho,ei,ej,nedge,n);

  M1(b1,b2,i)=gt11;
  M2(b1,b2,i)=gt22;
        end
       
       
    end
end
    %%
    thresh=n/3;
   M=zeros((fin1-strt1)/res1,(fin2-strt2)/res2);
    for b11=1:(fin1-strt1)/res1
        for b22=1:(fin2-strt2)/res2
            
            good=find(M1(b11,b22,:)+M2(b11,b22,:)>thresh);
            if length(good)>10
                M(b11,b22)=sum(M1(b11,b22,good))/sum(M2(b11,b22,good));
            end
        end
    end
    figure();
    hold on;
    [hwin1,hwin2]=find(M<1);
    hwins=[hwin1.*res1+strt1,hwin2.*res2+strt2];
    plot(hwins(:,1),hwins(:,2),'r*')
    [lwin1,lwin2]=find(M>1);
    lwins=[lwin1.*res1+strt1,lwin2.*res2+strt2];
    plot(lwins(:,1),lwins(:,2),'g*')
    [failb1,failb2]=find(M==0);
    fails=[failb1.*res1+strt1,failb2.*res2+strt2];
    plot(fails(:,1),fails(:,2),'k*')
