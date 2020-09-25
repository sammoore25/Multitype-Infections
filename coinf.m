N=1000;
c=4;

A=double(triu(rand(N)<c/N,1)); A=sparse(A+A');

edges=find(A); 
M=length(edges);
AE=sparse(N,N);
AE(edges)=(1:M);   
[I,J]=find(A);
    
H=sparse(M,M);
for e=1:M
    H(e,AE(J(e),setdiff(find(A(J(e),:)),I(e))))=1;
end

%%

F=[0.8,0.1,0.1;0.1,0.5,0.4;0.1,0.5,0.4];

Samples=100;


t=cell(3);
for a=1:3
    for b=1:3
        T=sparse(N,N);
        T(edges)=exprnd(1/F(a,b),M,1);
        t{a,b}=T;
    end
end

        
%%
        
        

for lambda=0.1:0.1:1
    
    r=1/lambda-1;
    
    tr=kron(ones(3,1),exprnd(1/r,N,1));
    
    
end
            
            
    



