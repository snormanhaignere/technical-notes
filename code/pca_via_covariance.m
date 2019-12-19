% Computation of PCAs from iterative computation of covariances

% dimesionality and data matrix
N = 3;
P = 4;
reps = 3;
X = randn(N*reps, P);

% PCA via svd
[U, S, V] = svd(X, 'econ');

% iterative computation of covariances
C = zeros(P, P);
for i = 1:reps
    xi = (1:N)+N*(i-1);
    C = C+X(xi, :)' * X(xi,:);
end

% principal component weights from eigenvectors of covariance
[pca_weights, pca_eigvals] = eig(C);

% sort by eigenvalues
pca_eigvals = diag(pca_eigvals);
[~, xi] = sort(pca_eigvals, 'descend');
pca_weights = pca_weights(:,xi);
pca_eigvals = pca_eigvals(xi);

% pca responses
pca_responses = X * pca_weights;

% compare
clc;
pca_responses
US = U*S

%% Linear regression using PCA to get the true weights

clc
wtrue = randn(P,1)
y = X*wtrue + randn(N*reps,1)*0;
west = pinv(X)*y
wpca = pinv(pca_responses(:,4))*y
w_pca_est = pca_weights(:,4) * wpca
