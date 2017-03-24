D = 10;
S = randn(D,1);
N1 = randn(D,1);
N2 = randn(D,1);

% X1, X2 are measurements
% S is a signal
% N1, N2 are iid samples of Gaussian noise
X1 = S + N1;
X2 = S + N2;

% variances of measurements
var_X1 = var(X1,1);
var_X2 = var(X2,1);
var_E = var(X1 - X2,1);

% estimates of the variances of the noise and signal
var_S = (var_X1 + var_X2 - var_E)/2;
% var_S = (var_X1*2 - var_E)/2;
var_N = var_E/2;

% correlation
C_X1_X2 = var_S / (var_S + var_N)
corr(X1,X2)