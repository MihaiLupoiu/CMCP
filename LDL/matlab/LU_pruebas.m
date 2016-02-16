clear;

A = [[2,2,1];[2,5,2];[1,2,2]];
%[L,D] = ldl(A);
%L';

n=3;
L=eye(n);
for k=1:n
    for j=k+1:n
        L(j,k)=A(j,k)/A(k,k);
    end
    for j=k+1:n
        for i=1:n
            A(j,i) = A(j,i)-L(j,k)*A(k,i);
        end
    end
end



%% Generar Matriz Simetrica
%B = triu(randn(100));
%A=B+B'-diag(diag(B));

