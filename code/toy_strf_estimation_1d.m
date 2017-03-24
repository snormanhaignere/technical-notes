clc;

% stimulus and impulse response
N = 5;
h = rand(N,1);
s = rand(N,1);

%% circular convolution by hand

r = zeros(N,1);
for t = 0:N-1
    for i = 0:N-1
        r(t+1) = r(t+1) + h(i+1)*s(mod(t-i,N)+1);
    end
end
r

%% circular convolution via fft

r = ifft(fft(h) .* fft(s))

%% convolution via matrix multiplication

S = nan(N,N);
for i = 1:N
    S(:,i) = circshift(s,i-1);
end
r = S * h

%% deconvolution via pseudoinverse

hh = pinv(S) * r
h

%% deconvolution via eigenvalue decomposition

clc;

% diagonalization
Css = S'*S
[Q,D] = eig(Css);
QDQ = Q * D * Q'

% solution via diagonalization
Csr = S' * r;
hh = Q * (D \ Q') * Csr
h


%% eigenvalue decomposition through DFT

clc;

% dft matrix
Qdft = dftmtx(N) / sqrt(N);

% stimulus power spectrum
% equal to fourier transform of stimulus autocorrelation function
Ps = abs(sqrt(N) * Qdft * s).^2
Fss = fft(S' * S(:,1))

% power spectrum also equals eigenvalues of the stimulus autocorrelation matrix
Css = S' * S;
Ddft = diag(Ps);
Css * Qdft
Qdft * Ddft

%% deconvolution through dft

% DFT thus gives alternative decomposition to that returned by eig
Css = S' * S;
QDQ = Qdft * Ddft * Qdft'

% solution via diagonalization
Csr = S' * r;
hh = Qdft * (Ddft \ Qdft') * Csr
h

%% deconvolution via power spectral decomposition

clc;

% fft of cross-correlation
srev = circshift(flipud(s), -(N-1));
Fsr = fft(srev) .* fft(r) % Dft of the cross-correlation function
fft(S' * r)

% power spectrum of stimulus
% equal to FFT of autocorrelation function
Fss = abs(fft(s)).^2;

% divide frequency spectrum
hh = ifft(Fsr ./ Fss)
h


%%

r = (1:5)'

% fft of cross-correlation
srev = circshift(flipud(s),1)
ifft(fft(srev) .* fft(r)) % Dft of the cross-correlation function
S' * r

rrev= circshift(flipud(r),1)
ifft(fft(s) .* fft(rrev))



