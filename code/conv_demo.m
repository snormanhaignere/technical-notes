t = 0:1:10;
result = nan(size(t));
close all;
figh = figure;
set(figh, 'Position', [100 100 400 800]);

for timepoint = [0:9]
    
    % signal
    signal = zeros(size(t));
    signal(t>2 & t<7) = 1;
    
    % flipped impulse response
    step_interval = [0,1; 1,2]+timepoint;
    h = zeros(size(t));
    h(t>=step_interval(1,1) & t<step_interval(1,2)) = -1;
    h(t>=step_interval(2,1) & t<step_interval(2,2)) = 1;
    
    % result
    result(t==(timepoint+1)) = sum(signal.*h);
    
    plot_style = {'k-o', 'LineWidth', 2};

    subplot(4,1,1);
    stem(t, signal, plot_style{:});
    ylim([-1 1]*2);
    xlabel('Time');
    title('Signal');
    
    subplot(4,1,2);
    stem(t, h, plot_style{:});
    ylim([-1 1]*2);
    title('Impulse Response Function (IRF)');
    xlabel('Time');
    
    subplot(4,1,3);
    stem(t, signal .* h, plot_style{:});
    ylim([-1 1]*2);
    title('Signal x IRF');
    xlabel('Time');
    
    subplot(4,1,4);
    stem(t, result, plot_style{:});
    ylim([-1 1]*2);
    xlim(t([1,end]))
    title('sum(Signal x IRF)');
    xlabel('Time');
    
    % get input for next time step
    input('Press enter for next time step');
    
end