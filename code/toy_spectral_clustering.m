x = [2,3,2,4,2,3];
n = length(x);
A = double(abs(x'*ones(1,n) - ones(n,1)*x)>0);
D = diag(sum(A,2));
L = sqrt(D)*A*sqrt(D);
[U,S,V] = svd(L);

[U,S,V] = svd(A);
