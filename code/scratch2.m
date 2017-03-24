% script to play with the FT of toy signals

addpath(genpath('/Users/svnh2/Dropbox (MIT)/mcdexp-svnh/general-analysis-code'));


T = 6;
F = 6;
% FTX = fft2(randn(T,F));
% FTX(1,1) = 0;
FTX = zeros(T,F);
FTX(1,4) = 1;
% FTX(:,6) = 1; 
% FTX(6,2) = 1;

% zero 2nd and fourth quadrant
t_freqs = fft_freqs_from_siglen(T, 1) * ones(1,F);
s_freqs = ones(T,1) * fft_freqs_from_siglen(F, 1)';


fftshift(t_freqs)
fftshift_nyqlast(s_freqs)

% % zero first and third, or second and fourth quadrant
% % depending on the sign of the temporal modulation rate
% FTX(sign(t_freqs) == -1 & s_freqs > 0) = 0;
% FTX(sign(t_freqs) == 1 & s_freqs < 0) = 0;

% zero first and third, or second and fourth quadrant
% depending on the sign of the temporal modulation rate
% FTX(~(sign(t_freqs) >= 0 & s_freqs >= 0)) = 0;


clc;
fftshift(t_freqs)
fftshift(s_freqs)
fftshift(abs(FTX))
% imagesc(fftshift(abs(FTX))


% plot magnitude
t = fftshift_nyqlast(t_freqs(:,1));
s = fftshift_nyqlast(s_freqs(1,:));

figure;
set(gcf, 'Position', [0 0 1000 200]);
subplot(1,3,1);
imagesc(s', t, fftshift_nyqlast(abs(FTX)));
set(gca,'YDir','normal', 'XTick', s, 'YTick', t);
xlabel('Freq (Pi Radians)'); ylabel('Freq (Pi Radians)');
colorbar;
title('Magnitude');
subplot(1,3,2);
imagesc(real(ifft2(FTX)))
subplot(1,3,3);
imagesc(imag(ifft2(FTX)))

%%
figure;
set(gcf, 'Position', [0 0 1000 200]);
subplot(1,3,1);
imagesc(s', t, fftshift_nyqlast(abs(fft2(real(ifft2(FTX))))));
set(gca,'YDir','normal', 'XTick', s, 'YTick', t);
xlabel('Freq (Pi Radians)'); ylabel('Freq (Pi Radians)');
colorbar;
title('Magnitude');
subplot(1,3,2);
imagesc(real(ifft2(FTX)))
subplot(1,3,3);
imagesc(imag(ifft2(FTX)))


%%

% FT1 = [0,randn(1,2) + sqrt(-1)*randn(1,2)];
% FT1 = [FT1, conj(flip(FT1(2:end)))]';
% FT2 = [0,randn(1,2) + sqrt(-1)*randn(1,2)];
% FT2 = [FT2, conj(flip(FT2(2:end)))]';
% 
% FTX = FT1 * FT2';

% ifft2(FT1 * FT2')


% FTX(1:3,2) = 1;
% FTX([1,4:5],5) = 1;
% f = fftshift(fft_freqs_from_siglen(N,1));
% if mod(N,2)==0
%     f(1) = -f(1);
% end
% 
% ifft2(FTX)


figure;
set(gcf, 'Position', [0 0 1000 200]);
subplot(1,3,1);
imagesc(abs(fftshift(FTX)));
set(gca,'YDir','normal');
xlabel('Freq (Pi Radians)'); ylabel('Freq (Pi Radians)');


%%
imagesc(fftshift(t_freqs));


