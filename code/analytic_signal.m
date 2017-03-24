f = [0 0.25 0.5 -0.25]
x = randn(4,1)';
FFT_x = fft(x)
FFT_hilbert = fft(hilbert(x))
