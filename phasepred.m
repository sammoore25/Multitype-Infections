%phase diagram - prediction


n=1000;


T=500;
c=3;
p=c/n;
rho=1;
mu1=1/10;
mu2=1/3;

its=1000;






% 
G = sparse(rand(n,n)) < p;
G = sparse(triu(G,1));
G = G + G';


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

lamb=eigs(hash,1);

eigt=@(bet1,bet2)teigs(bet1,bet2,rho,mu1,mu2);

finp1=0.6;resp1=0.002;strtp1=0;
finp2=2.2;resp2=0.005;strtp2=0;
MP=zeros((finp1-strtp1)/resp1,(finp2-strtp2)/resp2);
MP2=zeros((finp1-strtp1)/resp1,(finp2-strtp2)/resp2);
for b1=1:(finp1-strtp1)/resp1
for b2=1:(finp2-strtp2)/resp2
bet1=b1*resp1+strtp1;
bet2=b2*resp2+strtp2;
[e,E]=eigt(bet1,bet2);
MP(b1,b2)=e-1/lamb;
MP2(b1,b2)=E(1)-E(2);
end
end


%%
figure();hold on;
[hwin1,hwin2]=find(MP2<0);
hwins=[hwin1.*resp1+strtp1,hwin2.*resp2+strtp2];
plot(hwins(:,1),hwins(:,2),'.','color','[1 0.7 0.7]')
[lwin1,lwin2]=find(MP2>0);
lwins=[lwin1.*resp1+strtp1,lwin2.*resp2+strtp2];
plot(lwins(:,1),lwins(:,2),'.','color','[0.7 1 0.7]')
[failb1,failb2]=find(MP<0);
fails=[failb1.*resp1+strtp1,failb2.*resp2+strtp2];
plot(fails(:,1),fails(:,2),'.','color','[0.4 0.4 0.4]')
