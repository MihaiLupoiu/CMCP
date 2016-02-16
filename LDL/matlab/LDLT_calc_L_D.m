clear;
A = [[2,2,1];[2,5,2];[1,2,2]];
%[L,D] = ldl(A);

n=3;
D=zeros(n,1);
L=zeros(n,n);
v=zeros(n,1);

for j=1:n
    for i=1:j-1
        v(i) = D(i)*L(j,i);
    end
    
    ts=0;
    for k=1:j-1
        ts=ts+L(j,k)*v(k);
    end
    
    v(j)=A(j,j)-ts;
    D(j) = v(j);
    
    for i=j+1:n
        ts=0;
        for k=1:j-1
            ts=ts+L(i,k)*v(k);
        end
        L(i,j)=(A(i,j)-ts)/v(j);
    end    
end