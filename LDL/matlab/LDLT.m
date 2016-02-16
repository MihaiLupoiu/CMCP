clear;

A = [[2,2,1];[2,5,2];[1,2,2]];
%[L,D] = ldl(A);
%L';

n=3;
v=zeros(n,1);

for j=1:n
    for i=1:j-1
        v(i)=A(i,i)*A(j,i);
        v(i)
    end
    
    ts=0;
    for k=1:j-1
        ts=ts+A(j,k)*v(k);
    end
    
    v(j)=A(j,j)-ts;
    A(j,j)=v(j);
    
    for i=j+1:n
        ts=0;
        for k=1:j-1
            ts=ts+(A(i,k)*v(k));
        end
        A(i,j)=(A(i,j)-ts)/v(j);
    end
end



%% Generar Matriz Simetrica
%B = triu(randn(100));
%A=B+B'-diag(diag(B));

