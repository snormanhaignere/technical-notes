% This code illustrates how to compute variances, covariances and normalized
% covariances from the correlation and mean of the raw data

clc
N = 10;
X = rand(N, 3);

% mean and correlation/power of un-demeaned data
M = (1/N)*sum(X);
C = (1/N)*(X'*X);
P = diag(C)';

% variance from power/correlation and squared mean
X_demean = bsxfun(@minus, X, M);
V = (1/N) * sum(X_demean.^2)
V = var(X, 1)
V = P - M.^2

% covariance from power/correlation and mean outer-product
C_demean = (1/N) * X_demean' * X_demean
C_demean = cov(X, 1)
C_demean = C - M' * M

% code for z-scoring
S = sqrt(V);
X_zscore = bsxfun(@times, X_demean, 1./S)
X_zscore = zscore(X,1)

% code for normalized covariance
C_zscore = (1/N) * X_zscore' * X_zscore
C_zscore = 1./(S' * S) .* (C - M'*M)

%% The code below generalizes the above equations to complex numbers

clc
N = 10;
X = rand(N, 3) + rand(N,3)*sqrt(-1);

% mean and correlation/power of un-demeaned data
M = (1/N)*sum(X);
C = (1/N)*(X'*X);
P = diag(C)';

% % variance from power/correlation and squared mean
X_demean = bsxfun(@minus, X, M);
V = (1/N) * sum(conj(X_demean) .* X_demean)
V = var(X, 1)
V = P - conj(M) .* M

% covariance from power/correlation and mean outer-product
C_demean = (1/N) * X_demean' * X_demean
C_demean = cov(X, 1)
C_demean = C - M' * M

% code for z-scoring
S = sqrt(V);
X_zscore = bsxfun(@times, X_demean, 1./S)
X_zscore = zscore(X,1)

% code for normalized covariance
C_zscore = (1/N) * X_zscore' * X_zscore
C_zscore = 1./(S' * S) .* (C - M'*M)