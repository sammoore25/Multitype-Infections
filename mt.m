N=1000;
c=5;

delta=1;

A=double(triu(rand(N)<c/N,1)); A=sparse(A+A');

beta=0:0.005:0.2;
L=length(beta);

S=10;

n1=zeros(S,L);
n2=zeros(S,L);

for l=1:length(beta)
    
    for s=1:S
        
        T=[0.9,0.1;4,6]*beta(l);
        
        X=zeros(N,1);
        n=randi(N);
        X(n)=randi(2);
        
        I=1;
        inf=[n];
        inf2=[];
        
        while I>0
            
            inf2=[];
            
            for n=1:I
                
                i=inf(n);
                
                r=-log(rand)/delta;
                
                J=find(A(i,:));
                
                for m=1:length(J)
                    
                    j=J(m);
                    
                    if X(j)==0
                        
                        t=-log(rand(2,1))./T(X(i),:);
                        
                        if (t(1)<r)&(t(1)<t(2))
                            
                            inf2=[inf2,j];
                            X(j)=1;
                            
                        elseif (t(2)<r)&(t(2)<t(1))
                            
                            inf2=[inf2,j];
                            X(j)=2;
                            
                        end
                        
                    end
                    
                end
                
                X(i)=-X(i);
                
            end
            
            inf=inf2;
            
            I=length(inf);
            
        end
        
        n1(s,l)=sum(X==-1)/N;
        n2(s,l)=sum(X==-2)/N;
        
    end
end

%%
plot(beta,mean(n1)./mean(n1>0.02),'o',beta,mean(n2)./mean(n2>0.02),'o');



