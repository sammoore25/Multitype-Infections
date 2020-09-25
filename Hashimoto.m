function H=Hashimoto(A)

    N=length(A);
    edges=find(A); 
    M=length(edges);
    A(edges)=(1:M);   
    [I,J]=find(A);
    
    H=sparse(M,M);
    for e=1:M
        H(e,A(J(e),setdiff(find(A(J(e),:)),I(e))))=1;
    end
    
end
        
        
    
    
    
    