% Shows the relationship between ppca and svd.m
% 
% 2019-12-19: Created, Sam NH

clc;

% demean for equivalency
N = 6;
P = 3;
X = rand(N, P);
X = bsxfun(@minus, X, mean(X));

% run ppca and svd
K = P-1;
[Vp,Up] = ppca(X,K);
[U,S,V] = svd(X, 'econ');

% illustrate properties
fprintf('Up equals U*S:\n');
Up %#ok<NOPTS>
US = U*S %#ok<NOPTS>
fprintf('Vp equals V:\n');
Vp
V