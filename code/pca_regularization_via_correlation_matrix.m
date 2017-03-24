% Demonstrates that computing the correlation of a PCA-regularized matrix is equivalent to 
% reconstructing the correlation matrix directly from its top N PCs

P_total = 6;
P_to_keep = 3;
X = rand(10,P_total);
Cxx = X' * X

%% Regularize correlation matrix directly

[Q,S] = eig(Cxx);

% sort by eigen values
eigvals = diag(S);
[~,xi] = sort(eigvals,'descend');
Q = Q(:,xi);
S = diag(eigvals(xi));
eigvals = eigvals(xi);

% select top N% of eigen vectors
Q = Q(:,1:P_to_keep);
S = S(1:P_to_keep,1:P_to_keep);
Cxx_reg_direct = Q * S * Q'

%% Select top PCs of the original matrix and then compute correlation matrix

[U,S,V] = svd(X,'econ');
X_reg = U(:,1:P_to_keep) * S(1:P_to_keep,1:P_to_keep) * V(:,1:P_to_keep)';
Cxx_reg_via_pca = X_reg' * X_reg


