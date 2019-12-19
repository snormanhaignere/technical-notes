function conv_demo_gammatone(fc_Hz)

addpath(genpath('/Users/svnh2/Dropbox (MIT)/mindhive/spectrotemporal-synthesis-v2'));

% filter frequency
% try 150, 400, 6000
% fc_Hz = 6000;

% read in a snippet of speech
[wav,sr] = audioread('speech.wav');
win = [0,0.35]+0.15;
t = (0:length(wav)-1)/sr;
xi = t>=win(1) & t<win(2);
wav_snippet = wav(xi);
t = t(xi)-win(1);
wav_bounds = [min(wav_snippet), max(wav_snippet)];

% gammatone filter in frequency domain
N = 3*sr/fc_Hz;
H_filt = filt_temp_mod(fc_Hz, N, sr, false, false, true, 1, ...
    'gammatone', false, false, 1);

% analytic signal to get hilbert pair in real and imaginary
H_anal = analytic_from_spectrum(H_filt);

% convert to time domain
h_anal = ifft(H_anal);

% separate real and imaginary
h_real = real(h_anal);
h_imag = imag(h_anal);

% flip for convolution
h_real = flipud(h_real);
h_imag = flipud(h_imag);

scale_factor = max(abs(h_real));
h_real = h_real / scale_factor;
h_imag = h_imag / scale_factor;

figh = figure;
set(figh, 'Position', [100 100 600 800]);
result = zeros(size(wav_snippet));
env = zeros(size(wav_snippet));
started = false;
for timepoint = 1:(length(wav_snippet)-length(h_real))
    
    % place filter in time
    h_real_pad = zeros(size(wav_snippet));
    h_real_pad((1:length(h_real))+(timepoint-1)) = h_real;
    h_imag_pad = zeros(size(wav_snippet));
    h_imag_pad((1:length(h_imag))+(timepoint-1)) = h_imag;
    
    % calculate result of convolution with real
    result(timepoint + length(h_real)) = sum(wav_snippet .* h_real_pad);
    
    % calculate envelope using real and imaginary component
    env(timepoint + length(h_real)) = sqrt(sum(wav_snippet .* h_real_pad)^2+sum(wav_snippet .* h_imag_pad)^2);
    
    % plot every 100 steps
    
    if mod((timepoint-1),100)==0
        
        
        
        plot_style = {'k-', 'LineWidth', 1};
        
        % plot wav snippet
        subplot(4,1,1);
        plot(t, wav_snippet);
        xlabel('Time');
        ylim(wav_bounds);
        title('Signal');
        
        % plot real impulse response
        subplot(4,1,2);
        plot(t, h_real_pad, plot_style{:});
        ylim([-1,1]*max(abs(h_real_pad))*1.1);
        title('Impulse Response Function (IRF)');
        xlabel('Time');
        
        % plot the product of those two
        subplot(4,1,3);
        plot(t, wav_snippet .* h_real_pad, plot_style{:});
        ylim(wav_bounds);
        title('Signal x IRF');
        xlabel('Time');
        
        % plot the sum of the product along with the envelope
        subplot(4,1,4);
        plot(t, env, 'r-','LineWidth', 1);
        hold on;
        plot(t, result, plot_style{:});
        ylim(wav_bounds/scale_factor);
        title('sum(Signal x IRF)');
        xlabel('Time');
        drawnow;
        
        pause(0.01);
        
        if ~started
            input('Press enter to start demo');
            started = true;
        end
        
    end
end

