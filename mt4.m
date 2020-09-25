d=1;

mx=0.1;
my=0.33;

bx_sim=0.025:0.05:0.7;
by_sim=0.25:0.5:3;

Lx=length(bx_sim);
Ly=length(by_sim);

N=1000;
c=3;
A=double(triu(rand(N)<c/N,1)); A=sparse(A+A');

Samples=100;
FX=zeros(Samples,1);
FY=zeros(Samples,1);

Qx=[];
Qy=[];
Qz=[];

nchar = fprintf('%f',0);

for lx=1:Lx
for ly=1:Ly    
    
    R=diag([bx_sim(lx),by_sim(ly)])*[(1-mx),mx;my,(1-my)];  
    
    for sample=1:Samples
        
        fprintf(repmat('\b', 1, nchar));    
        nchar = fprintf('%f:%f  %f',lx/Lx,ly/Ly,sample/Samples);
            
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
                
                [ir,~,rr]=find(d*(x+y));
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
                       
            FX(sample)=X;
            FY(sample)=Y;         
            
    end
        
    if sum(FX)+sum(FY)<10*Samples
        Qz=[Qz;bx_sim(lx),by_sim(ly)];
    elseif sum(FX)>sum(FY)
        Qx=[Qx;bx_sim(lx),by_sim(ly)];
    else
        Qy=[Qy;bx_sim(lx),by_sim(ly)];
    end
end

    
end

%%

bx=0:0.01:0.7;
by=0:0.1:0.3;

qxy=(bx*d-2*bx*d*mx)./(d+2*bx*mx-2*bx*my-2*d*my);
qz=-d*(bx+d+bx*(mx-1)*c)./(d+d*(my-1)*c-bx*(c-1)*(1+c*(mx+my-1)));

area(bx,[max(qz,0);max(qxy-qz,0);3*ones(size(qz))]');
hold on;
plot(Qx(:,1),Qx(:,2),'ob',Qy(:,1),Qy(:,2),'or',Qz(:,1),Qz(:,2),'ok'); 
hold off;
axis([0,0.7,0,3])

