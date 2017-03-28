% Assesses the equations we typically use to estimate explained variance via
% regression

% external code repository
root_directory = my_root_directory;
addpath(genpath([root_directory '/general-analysis-code']));
addpath(genpath([root_directory '/export_fig_v2']));

% number of features
n_features = 5;

% number of times to repeat the analysis
n_reps = 100;
    
% factor to multiply the noise by
noise_factor = 1.5;

% metric used to quantify similarity
metrics = {'pearson-nofolds', 'pearson',  'demeaned-squared-error'};
n_metrics = length(metrics);

% parameter that determine the degree of nonlinearity
% higher values indicate greater nonlinearitiy
% alpha = linspace(1,3,10);
alpha = [1 1.25 1.5 1.75];
n_alpha = length(alpha);

% number of samples
sample_sizes = 50;
n_sample_sizes = length(sample_sizes);

% methods
% methods = {'least-squares'};
methods = {'ridge'};
n_methods = length(methods);

% number of folds to use for cross-validation
n_folds = 2;

% directory to save analysis results to
analysis_directory = [root_directory '/technical-notes/analysis' ...
    '/test-corrected-similarity-metrics'];
if ~exist(analysis_directory, 'dir'); mkdir(analysis_directory); end

% directory to save figures to
figure_directory = strrep(analysis_directory, 'analysis', 'figures');
if ~exist(figure_directory, 'dir'); mkdir(figure_directory); end

true_r = nan(n_reps, n_alpha, n_sample_sizes, n_methods, n_metrics);
noisy_r = nan(n_reps, n_alpha, n_sample_sizes, n_methods, n_metrics);
noise_corrected_r = nan(n_reps, n_alpha, n_sample_sizes, n_methods, n_metrics);
true_mse = nan(n_reps, n_alpha, n_sample_sizes, n_methods, n_metrics);
noisy_mse = nan(n_reps, n_alpha, n_sample_sizes, n_methods, n_metrics);
signal_test_retest = nan(n_reps, n_alpha, n_sample_sizes, n_methods, n_metrics);
prediction_test_retest = nan(n_reps, n_alpha, n_sample_sizes, n_methods, n_metrics);
for z = 1:length(metrics)
    
    switch metrics{z}
        case {'pearson', 'pearson-v1', 'pearson-nofolds'}
            simfunc = @fastcorr;
        case 'demeaned-squared-error'
            simfunc = @corr_variance_sensitive_symmetric;
        otherwise
            error('Switch statement fell through');
    end
    
    for q = 1:n_methods
        for k = 1:n_sample_sizes
            for i = 1:n_alpha
                
                % nonlinear function
                ResetRandStream2(1);
                b = randn(n_features,1);
                if alpha == 0
                    f_nonlin = @(f)(f * b);
                else
                    f_nonlin = @(f)(sin(f*alpha(i)) * b);
                end
                fprintf('%s, %s, sample size=%d, alpha=%d\n', ...
                    metrics{z}, methods{q}, sample_sizes(k), alpha(i)); drawnow;
                
                % string with the key parameter values
                param_idstring = ...
                    [metrics{z} '-' methods{q} ...
                    '-sampsize' num2str(sample_sizes(k)) ...
                    '-alpha' num2str(alpha(i)) ...
                    '-noisefactor' num2str(noise_factor) ...
                    '-nfeats' num2str(n_features) '-nfolds' num2str(n_folds) ...
                    '-nreps' num2str(n_reps)];
                
                MAT_file = [analysis_directory '/' param_idstring '.mat'];
                if ~exist(MAT_file, 'file')
                    S.true_r = nan(n_reps, 1);
                    S.true_mse = nan(n_reps, 1);
                    S.noisy_r = nan(n_reps, 1);
                    S.noisy_mse = nan(n_reps, 1);
                    S.noise_corrected_r = nan(n_reps, 1);
                    for j = 1:n_reps
                        
                        fprintf('Rep %d\n', j); drawnow;
                        
                        % sample features and compute function values
                        F = randn(sample_sizes(k),n_features);
                        y = f_nonlin(F);
                        
                        
                        % [w, bestK, mse, r] = regress_weights_from_2way_crossval(F(folds==2,:), y(folds==2), 2);
                        
                        % correlation with best linear approximation without noise
                        [linear_prediction, mse, r, folds] = regress_predictions_from_3way_crossval(...
                            F, y, n_folds, methods{q}, 2.^(-100:100), n_folds);
                        
                        % correlation for non-noisy data
                        S.true_r(j) = correlation_within_folds(...
                            y, linear_prediction, folds, metrics{z});
                        
                        % mean-squared error
                        S.true_mse(j) = mean((y-linear_prediction).^2);
                        
                        % two signals with distinct I.I.D. Gausian noise added
                        y_noisy1 = y + noise_factor * randn(sample_sizes(k),1);
                        y_noisy2 = y + noise_factor * randn(sample_sizes(k),1);
                        
                        % cross-validated predictions from noisy data
                        [noisy_prediction1, ~, ~, folds1] = regress_predictions_from_3way_crossval(...
                            F, y_noisy1, n_folds, methods{q}, 2.^(-100:100), n_folds);
                        [noisy_prediction2, ~, ~, folds2] = regress_predictions_from_3way_crossval(...
                            F, y_noisy2, n_folds, methods{q}, 2.^(-100:100), n_folds);
                        assert(all(folds1==folds2));
                        folds = folds1; clear folds1 folds2;
                        
                        % correlation of the noisy signal with the noisy prediction
                        S.noisy_r(j) = ...
                            simfunc(y_noisy1, noisy_prediction1)/2 ...
                            + simfunc(y_noisy2, noisy_prediction2)/2;
                        
                        % mean-squared error for noisy mse
                        S.noisy_mse(j) = mean(([y_noisy1;y_noisy2]...
                            -[noisy_prediction1;noisy_prediction2]).^2);
                        
                        % noise-correct
                        if strcmp(metrics{z}, 'pearson-v1')
                            S.noise_corrected_r(j) = ...
                                normalized_correlation_within_folds(...
                                [y_noisy1,y_noisy2], ...
                                [noisy_prediction1, noisy_prediction2], ...
                                folds);
                        elseif strcmp(metrics{z}, 'pearson-nofolds')
                            S.noise_corrected_r(j) = ...
                                normalized_correlation(...
                                [y_noisy1,y_noisy2], ...
                                [noisy_prediction1, noisy_prediction2]);
                        else
                            S.noise_corrected_r(j) = ...
                                normalized_correlation_within_folds_v2(...
                                [y_noisy1,y_noisy2], ...
                                [noisy_prediction1, noisy_prediction2], ...
                                folds, 'same_noise', false, 'metric', metrics{z});
                        end
                    end
                    
                    save(MAT_file, 'S');
                else
                    load(MAT_file, 'S');
                end
                
                true_r(:,i,k,q,z) = S.true_r;
                true_mse(:,i,k,q,z) = S.true_mse;
                noisy_r(:,i,k,q,z) = S.noisy_r;
                noisy_mse(:,i,k,q,z) = S.noisy_mse;
                noise_corrected_r(:,i,k,q,z) = S.noise_corrected_r;
            end
        end
    end
end

%% Plot correlation vs. alpha for fixed sample size

% n_reps, n_alpha, n_sample_sizes, n_methods, n_metrics

if n_alpha > 1
    
    for l = 1:n_sample_sizes
        figure;
        set(gcf, 'Position', [200 200 500*n_methods 500]);
        for q = 1:n_methods
            for z = 1:n_metrics
                
                subplot(n_metrics,n_methods, q + n_methods*(z-1));
                
                % -> rep x alpha x metric
                X = cat(ndims(true_r)+1, true_r, noisy_r, noise_corrected_r);
                X = squeeze_dims(X(:,:,l,q,z,:),[3,4,5]);
                
                % -> alpha x metric
                m = squeeze_dims(nanmedian(X,1),1);
                e = squeeze_dims(diff(quantile(X,[0.1587,1-0.1587],1)),1);
                
                % plot
                load('colormap-default-line-colors','cmap');
                h = nan(1,size(m,2));
                for i = 1:size(m,2);
                    h(i) = errorbar(alpha', m(:,i), e(:,i), ...
                        'o-', 'Color', cmap(i,:), 'LineWidth', 2);
                    hold on;
                end
                legend(h, {'clean', 'noisy', 'noise-corrected'}, 'Location', 'EastOutside');
                plot(xlim, [0 0], 'k--', 'LineWidth', 2);
                plot(xlim, [-1 -1], 'k--', 'LineWidth', 2);
                plot(xlim, [1 1], 'k--', 'LineWidth', 2);
                xlabel('Degree of Nonlinearity'); ylabel('Correlation (r)');
                title(sprintf('%s, %s', methods{q}, metrics{z}));
                ylim([-1.2 1.2]);
                bounds = [min(alpha(:)), max(alpha(:))];
                bounds = bounds + [-1 1]*diff(bounds)*0.1;
                xlim(bounds);
                fname = ['corr_vs_nonlinearity_noisefac' ...
                    num2str(noise_factor) '_smpsize' num2str(sample_sizes(l)) ...
                    '_nfeatures' num2str(n_features)];
                box off;
                set(gcf, 'PaperSize', [6*n_methods 6]);
                set(gcf, 'PaperPosition', [0.25 0.25 6*n_methods-0.5 5.5]);
                print([figure_directory '/' fname '.pdf'],'-dpdf')
                print([figure_directory '/' fname '.png'],'-dpng', '-r200');
            end
        end
    end
end

%% Plot mse vs. alpha for fixed sample size
% 
% if n_alpha > 1
%     
%     figure;
%     set(gcf, 'Position', [200 200 500*n_methods 500]);
%     for q = 1:n_methods
%         
%         subplot(1,n_methods,q);
%         
%         % -> rep x alpha x metric
%         X = log2(cat(5,true_mse, noisy_mse));
%         X = squeeze_dims(X(:,:,end,q,:),[3,4]);
%         
%         % -> alpha x metric
%         m = squeeze_dims(nanmedian(X,1),1);
%         e = squeeze_dims(diff(quantile(X,[0.1587,1-0.1587],1)),1);
%         
%         % plot
%         load('colormap-default-line-colors','cmap');
%         h = nan(1,size(m,2));
%         for i = 1:size(m,2);
%             h(i) = errorbar(alpha', m(:,i), e(:,i), 'Color', cmap(i,:), 'LineWidth', 2);
%             hold on;
%         end
%         legend(h, {'clean', 'noisy', 'noise-corrected'});
%         plot(xlim, [0 0], 'k--', 'LineWidth', 2);
%         xlabel('Degree of Nonlinearity'); ylabel('log2(Mean-Squared Error)');
%         fname = [pwd '/test_var_explained_mse_vs_nonlinearity_noiselevel' ...
%             num2str(noise_factor) '_smpsize' num2str(sample_sizes(end)) '_nfeatures' num2str(n_features)];
%         title(methods{q});
%         box off;
%         set(gcf, 'PaperSize', [6*n_methods 6]);
%         set(gcf, 'PaperPosition', [0.25 0.25 6*n_methods-0.5 5.5]);
%         print([fname '.pdf'],'-dpdf')
%         print([fname '.png'],'-dpng', '-r200');
%     end
% end

%% Plot correlation vs. sample size for fixed alpha

% n_reps, n_alpha, n_sample_sizes, n_methods, n_metrics

if n_sample_sizes > 1
    
    for l = 1:n_alpha
        figure;
        set(gcf, 'Position', [200 200 500*n_methods 500]);
        for q = 1:n_methods
            for z = 1:n_metrics
                
                subplot(n_metrics,n_methods, q + n_methods*(z-1));
                
                % -> rep x sample-size x metric
                X = cat(6, true_r, noisy_r, noise_corrected_r);
                X = squeeze_dims(X(:,l,:,q,z,:),[2 4 5]);
                
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
                legend(h, {'clean', 'noisy', 'noise-corrected'}, 'Location', 'EastOutside');
                plot(xlim, [0 0], 'k--', 'LineWidth', 2);
                plot(xlim, [1 1], 'k--', 'LineWidth', 2);
                plot(xlim, [-1 -1], 'k--', 'LineWidth', 2);
                xlabel('Sample Size'); ylabel('Correlation (r)');
                fname = ['corr_vs_smpsize_noisefac' ...
                    num2str(noise_factor) '_alpha' num2str(alpha(l)) ...
                    '_nfeatures' num2str(n_features)];
                title(sprintf('alpha %d, %s, %s', alpha(l), methods{q}, metrics{z}));
                box off;
                ylim([-1.2 1.2]);
                bounds = [min(log2(sample_sizes(:))), max(log2(sample_sizes(:)))];
                bounds = bounds + [-1 1]*diff(bounds)*0.1;
                xlim(bounds);
                set(gcf, 'PaperSize', [6*n_methods 6]);
                set(gcf, 'PaperPosition', [0.25 0.25 6*n_methods-0.5 5.5]);
                print([figure_directory '/' fname '.pdf'],'-dpdf')
                print([figure_directory '/' fname '.png'],'-dpng', '-r100');
            end
        end
    end
end

%% Plot mse vs. sample size for fixed alpha

% if n_sample_sizes > 1
%     
%     figure;
%     set(gcf, 'Position', [200 200 500*n_methods 500]);
%     for q = 1:n_methods
%         
%         subplot(1,n_methods,q);
%         
%         % -> rep x alpha x metric
%         X = log2(cat(5,true_mse, noisy_mse));
%         X = squeeze_dims(X(:,1,:,q,:),[2 4]);
%         
%         % -> alpha x metric
%         m = squeeze_dims(nanmean(X,1),1);
%         e = squeeze_dims(nanstd(X,1),1);
%         m = squeeze_dims(nanmedian(X,1),1);
%         e = squeeze_dims(diff(quantile(X,[0.1587,1-0.1587],1)),1);
%         
%         % plot
%         load('colormap-default-line-colors','cmap');
%         h = nan(1,size(m,2));
%         for i = 1:size(m,2);
%             h(i) = errorbar(log2(sample_sizes'), m(:,i), e(:,i), 'Color', cmap(i,:), 'LineWidth', 2);
%             hold on;
%         end
%         set(gca, 'XTick', log2(sample_sizes), 'XTickLabel', sample_sizes);
%         legend(h, {'clean', 'noisy'},'Location', 'Best');
%         plot(xlim, [0 0], 'k--', 'LineWidth', 2);
%         xlabel('Sample Size'); ylabel('log2(Mean-Squared Error)');
%         fname = [pwd '/test_var_explained_mse_vs_smpsize_noiselevel' ...
%             num2str(noise_factor) '_alpha' num2str(alpha(1)) '_nfeatures' num2str(n_features)];
%         title(methods{q});
%         box off;
%         ylim([-10 10]);
%         %         ylim([-1 1]);
%         set(gcf, 'PaperSize', [6*n_methods 6]);
%         set(gcf, 'PaperPosition', [0.25 0.25 6*n_methods-0.5 5.5]);
%         print([fname '.pdf'],'-dpdf')
%         print([fname '.png'],'-dpng', '-r100');
%     end
% end
