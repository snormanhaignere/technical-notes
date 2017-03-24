
addpath(genpath('/Users/svnh2/Dropbox (MIT)/mcdexp-svnh/general-analysis-code'));


N = 6;
X = randn(N,N);

%% Shifted FFT

f = fftshift_nyqlast(fft_freqs_from_siglen(N,1));
FTX = fftshift_nyqlast(fft2(X));

figure;
set(gcf, 'Position', [0 0 1000 300])

% plot magnitude
subplot(1,2,1);
imagesc(f, f, abs(FTX));
set(gca,'YDir','normal', 'XTick', f, 'YTick', f);
xlabel('Freq (Pi Radians)'); ylabel('Freq (Pi Radians)');
colorbar;
title('Magnitude');

% plot phase
subplot(1,2,2);
imagesc(f, f, angle(FTX));
set(gca,'YDir','normal', 'XTick', f, 'YTick', f);
xlabel('Freq (Pi Radians)'); ylabel('Freq (Pi Radians)');
colorbar;
title('Phase');

% color map
colormap(cbrewer('div', 'RdBu', 64, 'pchip'));

% save
export_fig(['FFT-2D-' num2str(N) 'x' num2str(N) '.png'], ...
    '-png','-transparent', '-r100');

%% Unshifted

f = fft_freqs_from_siglen(N,1);
FTX = fft2(X);

figure;
set(gcf, 'Position', [0 0 1000 300])

% plot magnitude
subplot(1,2,1);
x = 1:length(f);
imagesc(x, x, abs(FTX));
set(gca, 'XTick', x, 'XTickLabel', f, 'YTick', x, 'YTickLabel', f);
xlabel('Freq (Pi Radians)'); ylabel('Freq (Pi Radians)');
colorbar;
title('Magnitude');

% plot phase
subplot(1,2,2);
imagesc(x, x, angle(FTX));
set(gca, 'XTick', x, 'XTickLabel', f, 'YTick', x, 'YTickLabel', f);
xlabel('Freq (Pi Radians)'); ylabel('Freq (Pi Radians)');
colorbar;
title('Phase');

% color map
colormap(cbrewer('div', 'RdBu', 64, 'pchip'));

% save
export_fig(['FFT-2D-' num2str(N) 'x' num2str(N) '-Unshifted.png'], ...
    '-png','-transparent', '-r100');

