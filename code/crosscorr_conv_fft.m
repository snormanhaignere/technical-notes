% stimulus and impulse response
N = 5;
s = rand(N,1);
h = rand(N,1)

S = nan(N,N);
for i = 1:N
    S(:,i) = circshift(s,i-1);
end

% convolution
S * h
r = ifft(fft(s) .* fft(h))

% cross-correlation
S' * h
ifft(conj(fft(s)) .* fft(h))