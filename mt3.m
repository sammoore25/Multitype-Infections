delta=2;

mx=0.01;
my=0.25;

phi=2:0.01:8;


L=length(phi);
w=zeros(2,L);
rho=zeros(1,L);

for l=1:L   
    R=[(1-mx),mx;phi(l)*my,phi(l)*(1-my)];  
    T=diag(1./(sum(R')+delta))*R;
    [u,rho(l)]=eigs(T',1);
    w(:,l)=u/sum(u);
end
figure();
plot(phi,w,'-');

%%

phi_emp=2:0.5:8;
L=length(phi_emp);

% N=1000;
% c=4;
% A=double(triu(rand(N)<c/N,1)); A=sparse(A+A');

Samples=100;
FX=zeros(Samples,L);
FY=zeros(Samples,L);

nchar = fprintf('%f',0);

for l=1:L
    
    R=[(1-mx),mx;phi_emp(l)*my,phi_emp(l)*(1-my)];  
    
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
            
            if X+Y>10
                success=1;
            end
           
            
        end
            FX(sample,l)=X;
            FY(sample,l)=Y;
    end
    
end

%%
%plot(beta,FX','ok',beta,FY','or');
plot(phi,w(1,:),'b',phi,w(2,:),'r',phi_emp,mean(FX./(FX+FY)),'ob',phi_emp,mean(FY./(FX+FY)),'or')
% 
% %%
% phi=1:0.01:4;
% L=length(phi);
% h=zeros(1,L);
% for l=1:L
%     R=[(1-mx),mx;phi(l)*my,phi(l)*(1-my)];  
%     T=diag(1./(sum(R')+delta))*R;
%     [u,lambda]=eigs(T',1);
%     w(:,l)=u/sum(u);
%     h(l)=1-lambda-lambertw(-c*lambda*exp(-c*lambda))/c;
% end
% subplot(2,1,1);
% plot(phi,(1-h).*w(1,:),phi,(1-h).*w(2,:),phi_emp,mean(FX)/N/2,'ob',phi_emp,mean(FY)/N/2,'or');
% subplot(2,1,2);
% plot(phi,(1-h),'k',phi_emp,mean(FX+FY)/N/2,'ok');
%         
%     



