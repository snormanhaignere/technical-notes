% 

clc;
ResetRandStream2(1);

M = 3;
N = 6;
X1 = randn(M,N);
X2 = randn(M,N);

% whiten
[U1,S1,V1] = svd(X1,'econ');
X1_white = U1*V1';
[U2,S2,V2] = svd(X2,'econ');
X2_white = U2*V2';

% verify whitening
% [U1_white, S1_white, V1_white] = svd(X1_white, 'econ')

% concatenate
Xcat = [X1_white,X2_white];
Xcat = [U1, U2]
% Xcat = X1_white/2 + X2_white/2;

[U, S, V] = svd(Xcat, 'econ');

S



