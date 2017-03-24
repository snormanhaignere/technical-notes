x = randn(3,1);
y = randn(3,1);
f = x * y';
fft2(x * y')
fft(x,[],1) * fft(y',[],2)


% imagesc(imag(fft2(f)))