% Illustrates complex PCA
% 
% 2017-05-24: Created, Sam NH

% complex stimulus
X = randn(4,5) + sqrt(-1)*randn(4,5);

% covariance
C = X*transpose(conj(X));

% eigen values of covariance matrix
[E, D] = eig(C);

% are the same as those computed by the svd up to a phase shift
[U, S, V] = svd(X, 'econ');
mag = abs(E' * U)
ph = angle(E' * U)

% reconstruct
Xh = U*S*V';
abs(X - Xh)

% projection onto U
SV = U' * X;
SV - S*V'

% projection onto V
US = X * V;
US - U*S;