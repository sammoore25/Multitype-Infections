
delta=1;

mx=0.01;
my=0.25;
phi=4;

%F=[0.99,0.01;1,3];
F=[(1-mx),mx;phi*my,phi*(1-my)];

beta=[0.001,[0.1:0.1:2]];
L=length(beta);
w=zeros(2,L);
rho=zeros(1,L);

for l=1:L
    R=F*beta(l);    
    T=diag(1./(sum(R')+delta))*R;
    [u,rho(l)]=eigs(T',1);
    w(:,l)=u/sum(u);
end

plot(beta,w,'-');

%%

N=1000;
c=2;
A=double(triu(rand(N)<c/N,1)); A=sparse(A+A');

Samples=100;
FX=zeros(Samples,L);
FY=zeros(Samples,L);

nchar = fprintf('%f',0);

for l=1:L
    
    R=F*beta(l);
    
    for sample=1:Samples
        
        fprintf(repmat('\b', 1, nchar));    
        nchar = fprintf('%f:%f',l/L,sample/Samples);
        
        success=0;
        while success==0
            
            n=randi(N);
            s=ones(N,1); s(n)=0;
            x=sparse(N,1); 
            x(n)=round(rand);
            %if rand<w(1,l)
               %x(n)=1;
            %end
            y=sparse(N,1); y(n)=1-x(n);
            
            I=1;
            X=0;
            Y=0;
            
%             s=ones(N,1);
%             x=sparse(N,1);
%             y=sparse(N,1);
%             n=randsample(N,2);
%             x(n(1))=1;
%             y(n(2))=1;
%             s(n)=0;
%             I=2;
%             X=0;
%             Y=0;
            
            while (I>0)
                
                [ir,~,rr]=find(delta*(x+y));
                [ix,~,rx]=find(s.*(R(1,1)*A*x+R(2,1)*A*y));
                [iy,~,ry]=find(s.*(R(1,2)*A*x+R(2,2)*A*y));
                
                rates=[sum(rr),sum(rx),sum(ry)];
                
                switch randsample(3,1,true,rates);
                    
                    case 1
                        
                        i=ir(randsample(length(rr),1,true,rr));
                        
                        X=X+x(i);
                        Y=Y+y(i);
                        x(i)=0;
                        y(i)=0;
                        I=I-1;
                        
                    case 2
                        
                        i=ix(randsample(length(rx),1,true,rx));
                        
                        x(i)=1;
                        s(i)=0;
                        I=I+1;
                        
                    case 3
                        
                        i=iy(randsample(length(ry),1,true,ry));
                        
                        y(i)=1;
                        s(i)=0;
                        I=I+1;
                        
                end
                
            end
            
            if X+Y>0
                success=1;
            end
           
            
        end
            FX(sample,l)=X;
            FY(sample,l)=Y;
    end
    
end

%%
%plot(beta,FX','ok',beta,FY','or');
plot(beta,mean(FX./(FX+FY)),'ob',beta,w(1,:),beta,mean(FY./(FX+FY)),'or',beta,w(2,:))

