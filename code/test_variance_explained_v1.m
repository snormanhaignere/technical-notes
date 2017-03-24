% Assesses the equations we typically use to estimate explained variance via
% regression

% external code repository
root_directory = '/mindhive/nklab/u/svnh';
addpath(genpath([root_directory '/general-analysis-code']));
addpath(genpath([root_directory '/export_fig_v2']));

% number of features
n_features = 5;

% number of samples
% sample_sizes = 1000;
sample_sizes = [10 25 50 100 250 500 1000];
n_sample_sizes = length(sample_sizes);

% parameter that determine the degree of nonlinearity
% higher values indicate greater nonlinearitiy
% alpha = linspace(1,3,10);
alpha = 1;
n_alpha = length(alpha);

% number of times to repeat the analysis
n_reps = 100;

% factor to multiply the noise by
noise_factor = 1;

% methods
% methods = {'least-squares'};
methods = {'least-squares','ridge'};
n_methods = length(methods);

% number of folds to use for cross-validation
n_folds = 10;

% MAT_file = ['test_var_explained_' DataHash({sample_sizes, alpha, methods, n_reps, noise_factor}) '.mat'];
MAT_file = ['test_var_explained_' DataHash({sample_sizes, alpha, methods, ...
    n_reps, noise_factor, n_features, n_folds}) '.mat'];
if ~exist(MAT_file, 'file');
    warning('OFF'); %#ok<WNOFF>
    true_r = nan(n_reps, n_alpha, n_sample_sizes, n_methods);
    noisy_r = nan(n_reps, n_alpha, n_sample_sizes, n_methods);
    noise_corrected_r = nan(n_reps, n_alpha, n_sample_sizes, n_methods);
    true_mse = nan(n_reps, n_alpha, n_sample_sizes, n_methods);
    noisy_mse = nan(n_reps, n_alpha, n_sample_sizes, n_methods);
    signal_test_retest = nan(n_reps, n_alpha, n_sample_sizes, n_methods);
    prediction_test_retest = nan(n_reps, n_alpha, n_sample_sizes, n_methods);
    for q = 1:n_methods
        for k = 1:n_sample_sizes
            for i = 1:n_alpha
                
                % nonlinear function
                ResetRandStream2(alpha(i));
                b = randn(n_features,1);
                if alpha == 0
                    f_nonlin = @(f)(f * b);
                else
                    f_nonlin = @(f)(sin(f*alpha(i)) * b);
                end
                fprintf('%s, sample size=%d, alpha=%d\n', ...
                    methods{q}, sample_sizes(k), alpha(i)); drawnow;
                
                for j = 1:n_reps
                    
                    % sample features and compute function values
                    F = randn(sample_sizes(k),n_features);
                    y = f_nonlin(F);
                    
                    % correlation with best linear approximation without noise
                    linear_prediction = regress_predictions_from_3way_crossval(...
                        F, y, n_folds, methods{q});
                    true_r(j,i,k,q) = corr(y, linear_prediction);
                    
                    % mean-squared error
                    true_mse(j,i,k,q) = mean((y-linear_prediction).^2);
                    
                    % two signals with distinct I.I.D. Gausian noise added
                    y_noisy1 = y + noise_factor * randn(sample_sizes(k),1);
                    y_noisy2 = y + noise_factor * randn(sample_sizes(k),1);
                    
                    % cross-validated predictions from noisy data
                    noisy_prediction1 = regress_predictions_from_3way_crossval(...
                        F, y_noisy1, n_folds, methods{q});
                    noisy_prediction2 = regress_predictions_from_3way_crossval(...
                        F, y_noisy2, n_folds, methods{q});
                    
                    % correlation of the noisy signal with the noisy prediction
                    noisy_r(j,i,k,q) = tanh(atanh(corr(y_noisy1, noisy_prediction1))/2 ...
                        + atanh(corr(y_noisy2, noisy_prediction2))/2);
                    
                    % mean-squared error for noisy mse
                    noisy_mse(j,i,k,q) = mean(([y_noisy1;y_noisy2]...
                        -[noisy_prediction1;noisy_prediction2]).^2);
                    
                    % test retest of signal and predictions
                    signal_test_retest(j,i,k,q) = corr(y_noisy1,y_noisy2);
                    prediction_test_retest(j,i,k,q) = corr(noisy_prediction1,noisy_prediction2);
                    
                    % noise-correct
                    if signal_test_retest(j,i,k,q) > 0 && prediction_test_retest(j,i,k,q) > 0
                        %                         noise_corrected_r(j,i) = noisy_r(j,i) ...
                        %                             / sqrt(y_test_retest * y_prediction);
                        noise_corrected_r(j,i,k,q) = normalized_correlation([y_noisy1,y_noisy2],...
                            [noisy_prediction1, noisy_prediction2]);
                        
                    end
                end
            end
        end
    end
    warning('ON'); %#ok<WNON>
    save(MAT_file)
else
    load(MAT_file);
end

%% Plot correlation vs. alpha for fixed sample size

if n_alpha > 1
    
    figure;
    set(gcf, 'Position', [200 200 500*n_methods 500]);
    for q = 1:n_methods
        
        subplot(1,n_methods,q);
        
        % -> rep x alpha x metric
        X = cat(5,true_r, noisy_r, noise_corrected_r);
        X = squeeze_dims(X(:,:,end,q,:),[3,4]);
        
        % -> alpha x metric
        m = squeeze_dims(nanmedian(X,1),1);
        e = squeeze_dims(diff(quantile(X,[0.1587,1-0.1587],1)),1);
        
        % plot
        load('colormap-default-line-colors','cmap');
        h = nan(1,size(m,2));
        for i = 1:size(m,2);
            h(i) = errorbar(alpha', m(:,i), e(:,i), 'Color', cmap(i,:), 'LineWidth', 2);
            hold on;
        end
        legend(h, {'clean', 'noisy', 'noise-corrected'});
        plot(xlim, [0 0], 'k--', 'LineWidth', 2);
        xlabel('Degree of Nonlinearity'); ylabel('Correlation (r)');
        fname = [pwd '/test_var_explained_corr_vs_nonlinearity_noiselevel' ...
            num2str(noise_factor) '_smpsize' num2str(sample_sizes(end)) '_nfeatures' num2str(n_features)];
        title(methods{q});
        box off;
        set(gcf, 'PaperSize', [6*n_methods 6]);
        set(gcf, 'PaperPosition', [0.25 0.25 6*n_methods-0.5 5.5]);
        print([fname '.pdf'],'-dpdf')
        print([fname '.png'],'-dpng', '-r200');
    end
end

%% Plot mse vs. alpha for fixed sample size

if n_alpha > 1
    
    figure;
    set(gcf, 'Position', [200 200 500*n_methods 500]);
    for q = 1:n_methods
        
        subplot(1,n_methods,q);
        
        % -> rep x alpha x metric
        X = log2(cat(5,true_mse, noisy_mse));
        X = squeeze_dims(X(:,:,end,q,:),[3,4]);
        
        % -> alpha x metric
        m = squeeze_dims(nanmedian(X,1),1);
        e = squeeze_dims(diff(quantile(X,[0.1587,1-0.1587],1)),1);
        
        % plot
        load('colormap-default-line-colors','cmap');
        h = nan(1,size(m,2));
        for i = 1:size(m,2);
            h(i) = errorbar(alpha', m(:,i), e(:,i), 'Color', cmap(i,:), 'LineWidth', 2);
            hold on;
        end
        legend(h, {'clean', 'noisy', 'noise-corrected'});
        plot(xlim, [0 0], 'k--', 'LineWidth', 2);
        xlabel('Degree of Nonlinearity'); ylabel('log2(Mean-Squared Error)');
        fname = [pwd '/test_var_explained_mse_vs_nonlinearity_noiselevel' ...
            num2str(noise_factor) '_smpsize' num2str(sample_sizes(end)) '_nfeatures' num2str(n_features)];
        title(methods{q});
        box off;
        set(gcf, 'PaperSize', [6*n_methods 6]);
        set(gcf, 'PaperPosition', [0.25 0.25 6*n_methods-0.5 5.5]);
        print([fname '.pdf'],'-dpdf')
        print([fname '.png'],'-dpng', '-r200');
    end
end

%% Plot correlation vs. sample size for fixed alpha

if n_sample_sizes > 1
    
    figure;
    set(gcf, 'Position', [200 200 500*n_methods 500]);
    for q = 1:n_methods
        
        subplot(1,n_methods,q);
        
        % -> rep x alpha x metric
        X = cat(5,true_r, noisy_r, noise_corrected_r);
        X = squeeze_dims(X(:,1,:,q,:),[2 4]);
                
        % -> alpha x metric
        m = squeeze_dims(nanmedian(X,1),1);
        e = squeeze_dims(diff(quantile(X,[0.1587,1-0.1587],1)),1);
                
        % plot
        load('colormap-default-line-colors','cmap');
        h = nan(1,size(m,2));
        for i = 1:size(m,2);
            h(i) = errorbar(log2(sample_sizes'), m(:,i), e(:,i), 'Color', cmap(i,:), 'LineWidth', 2);
            hold on;
        end
        set(gca, 'XTick', log2(sample_sizes), 'XTickLabel', sample_sizes);
        legend(h, {'clean', 'noisy', 'noise-corrected'},'Location', 'Best');
        plot(xlim, [0 0], 'k--', 'LineWidth', 2);
        xlabel('Sample Size'); ylabel('Correlation (r)');
        fname = [pwd '/test_var_explained_corr_vs_smpsize_noiselevel' ...
            num2str(noise_factor) '_alpha' num2str(alpha(1)) '_nfeatures' num2str(n_features)];
        title(methods{q});
        box off;
        ylim([-1 1]);
        set(gcf, 'PaperSize', [6*n_methods 6]);
        set(gcf, 'PaperPosition', [0.25 0.25 6*n_methods-0.5 5.5]);
        print([fname '.pdf'],'-dpdf')
        print([fname '.png'],'-dpng', '-r100');
    end
end

%% Plot mse vs. sample size for fixed alpha

if n_sample_sizes > 1
    
    figure;
    set(gcf, 'Position', [200 200 500*n_methods 500]);
    for q = 1:n_methods
        
        subplot(1,n_methods,q);
        
        % -> rep x alpha x metric
        X = log2(cat(5,true_mse, noisy_mse));
        X = squeeze_dims(X(:,1,:,q,:),[2 4]);
        
        % -> alpha x metric
        m = squeeze_dims(nanmean(X,1),1);
        e = squeeze_dims(nanstd(X,1),1);
        m = squeeze_dims(nanmedian(X,1),1);
        e = squeeze_dims(diff(quantile(X,[0.1587,1-0.1587],1)),1);
        
        % plot
        load('colormap-default-line-colors','cmap');
        h = nan(1,size(m,2));
        for i = 1:size(m,2);
            h(i) = errorbar(log2(sample_sizes'), m(:,i), e(:,i), 'Color', cmap(i,:), 'LineWidth', 2);
            hold on;
        end
        set(gca, 'XTick', log2(sample_sizes), 'XTickLabel', sample_sizes);
        legend(h, {'clean', 'noisy'},'Location', 'Best');
        plot(xlim, [0 0], 'k--', 'LineWidth', 2);
        xlabel('Sample Size'); ylabel('log2(Mean-Squared Error)');
        fname = [pwd '/test_var_explained_mse_vs_smpsize_noiselevel' ...
            num2str(noise_factor) '_alpha' num2str(alpha(1)) '_nfeatures' num2str(n_features)];
        title(methods{q});
        box off;
        ylim([-10 10]);
%         ylim([-1 1]);
        set(gcf, 'PaperSize', [6*n_methods 6]);
        set(gcf, 'PaperPosition', [0.25 0.25 6*n_methods-0.5 5.5]);
        print([fname '.pdf'],'-dpdf')
        print([fname '.png'],'-dpng', '-r100');
    end
end
