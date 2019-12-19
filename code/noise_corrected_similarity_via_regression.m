function noise_corrected_similarity_via_regression(...
    correction_method, average_before_combining_terms, variance_centering, ...
    n_folds, n_features, noise_factor, noise_corrected_crossval, regularization_metric) 

% Assesses the equations we typically use to estimate explained variance via
% regression
% 
% 2017-04-01: Several updates, SNH

% external code repository
root_directory = my_root_directory;
addpath(genpath([root_directory '/general-analysis-code']));
addpath(genpath([root_directory '/export_fig_v2']));

% number of features
% n_features = 10;

% number of times to repeat the analysis
n_reps = 100;
    
% factor to multiply the noise by
% noise_factor = 1.5;

% how to average across folds
% average_before_combining_terms = true;

% correction method
% correction_method = 'variance-based';

% whether or not to normalize variance across repetitions
% only relevant for variance-based correction
% variance_centering = false;

% metric used to quantify similarity
if strcmp(correction_method, 'correlation-based')
    metrics = {'pearson'};
else
    metrics = {'normalized-squared-error'};%, 'pearson', 'demeaned-squared-error'};
end
n_metrics = length(metrics);

% metric used for chosing the regularization parameter
% regularization_metric = 'unnormalized-squared-error';

% whether or not to use the noise-corrected metric to select the regularization 
% parameter
% noise_corrected_crossval = true;

% parameter that determine the degree of nonlinearity
% higher values indicate greater nonlinearitiy
% alpha = linspace(1,3,10);
alpha = [1 1.25 1.5 1.75 2];
n_alpha = length(alpha);

% number of samples
sample_sizes = 165;
n_sample_sizes = length(sample_sizes);

% methods
% methods = {'least-squares'};
methods = {'ridge'};
n_methods = length(methods);

% number of folds to use for cross-validation
% n_folds = 2;

% directory to save analysis results to
analysis_directory = [root_directory '/technical-notes/analysis' ...
    '/noise-corrected-similarity-via-regression'];
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
        case {'pearson'}
            funcname = 'pearson';
            simfunc = @fastcorr;
        case {'demeaned-squared-error'}
            funcname = 'demeaned-squared-error';
            simfunc = @corr_variance_sensitive_symmetric;
        case {'normalized-squared-error'}
            funcname = 'normalized-squared-error';
            simfunc = @normalized_squared_error;
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
                    '-nreps' num2str(n_reps) ...
                    '-averagebeforecombining' num2str(average_before_combining_terms) ...
                    '-varcenter' num2str(variance_centering) ...
                    '-' correction_method '-' regularization_metric ...
                    '-noisecorrCV' num2str(noise_corrected_crossval)];
                                
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
                        
                        % correlation with best linear approximation without noise
                        [linear_prediction, ~, ~, folds] = ...
                            regress_predictions_from_3way_crossval(...
                            F, y, 'test_folds', n_folds, 'train_folds', n_folds, ...
                            'method', methods{q}, 'K', 2.^(-100:100), ...
                            'regularization_metric', regularization_metric, ...
                            'warning', false);
                        
                        % correlation for non-noisy data
                        S.true_r(j) = correlation_within_folds(...
                            y, linear_prediction, folds, funcname);
                        
                        % mean-squared error
                        S.true_mse(j) = mean((y-linear_prediction).^2);
                        
                        % two signals with distinct I.I.D. Gausian noise added
                        y_noisy1 = y + noise_factor * randn(sample_sizes(k),1);
                        y_noisy2 = y + noise_factor * randn(sample_sizes(k),1);
                        
                        % cross-validated predictions from noisy data
                        if noise_corrected_crossval
                            [noisy_predictions, folds] = ...
                                regress_predictions_from_3way_crossval_noisecorr(...
                                F, [y_noisy1, y_noisy2], ...
                                'test_folds', n_folds, 'train_folds', n_folds, ...
                                'method', methods{q}, 'K', 2.^(-100:100), ...
                                'regularization_metric', regularization_metric, ...
                                'correction_method', correction_method, ...
                                'warning', false);
                            noisy_prediction1 = noisy_predictions(:,1);
                            noisy_prediction2 = noisy_predictions(:,2);
                        else
                            [noisy_prediction1, ~, ~, folds1] = ...
                                regress_predictions_from_3way_crossval(...
                                F, y_noisy1, ...
                                'test_folds', n_folds, 'train_folds', n_folds, ...
                                'method', methods{q}, 'K', 2.^(-100:100), ...
                                'regularization_metric', regularization_metric, ...
                                'warning', false);
                            [noisy_prediction2, ~, ~, folds2] = regress_predictions_from_3way_crossval(...
                                F, y_noisy2, 'test_folds', n_folds, ...
                                'test_folds', n_folds, 'train_folds', n_folds, ...
                                'method', methods{q}, 'K', 2.^(-100:100), ...
                                'regularization_metric', regularization_metric, ...
                                'warning', false);
                            assert(all(folds1==folds2));
                            folds = folds1; clear folds1 folds2;
                        end

                        % correlation of the noisy signal with the noisy prediction
                        S.noisy_r(j) = ...
                            correlation_within_folds(...
                            y_noisy1, noisy_prediction1, folds, funcname)/2 ...
                            + correlation_within_folds(...
                            y_noisy2, noisy_prediction2, folds, funcname)/2;
                        
                        % mean-squared error for noisy mse
                        S.noisy_mse(j) = mean(([y_noisy1; y_noisy2]...
                            -[noisy_prediction1; noisy_prediction2]).^2);
                        
                        switch correction_method
                            case 'variance-based'
                                S.noise_corrected_r(j) = ...
                                    noise_corrected_similarity_within_folds(...
                                    [y_noisy1,y_noisy2], ...
                                    [noisy_prediction1, noisy_prediction2], ...
                                    folds, 'same_noise', false, 'metric', metrics{z}, ...
                                    'variance_centering', variance_centering, ...
                                    'average_before_combining_terms', ...
                                    average_before_combining_terms);
                                
                            case 'correlation-based'
                                assert(strcmp(metrics{z}, 'pearson'));
                                S.noise_corrected_r(j) = ...
                                    normalized_correlation_within_folds(...
                                    [y_noisy1,y_noisy2], ...
                                    [noisy_prediction1, noisy_prediction2], ...
                                    folds, 'z_averaging', false, ...
                                    'average_before_combining_terms', ...
                                    average_before_combining_terms);
                                
                            otherwise
                                error('Switch statement fell through');
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
        for q = 1:n_methods
            for z = 1:n_metrics
                
                N = n_metrics*n_methods;
                n_rows = round(sqrt(N));
                n_cols = ceil(N/n_rows);
                set(gcf, 'Position', [200 200 400*n_cols 300*n_rows]);
                set(gcf, 'Color', [1 1 1]);
                subplot(n_rows,n_cols, q + n_methods*(z-1));
                
                % -> rep x alpha x metric
                X = cat(6, true_r, noisy_r, noise_corrected_r);
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

                %                 set(gcf, 'PaperSize', [6*n_methods 6]);
                %                 set(gcf, 'PaperPosition', [0.25 0.25 6*n_methods-0.5 5.5]);
                %                 print([figure_directory '/' fname '.pdf'],'-dpdf')
                %                 print([figure_directory '/' fname '.png'],'-dpng', '-r200');
            end
        end
        
        fname = ['corr-vs-nonlinearity-noisefac' ...
            num2str(noise_factor) '-smpsize' num2str(sample_sizes(l)) ...
            '-nfeatures' num2str(n_features) ...
            '-averagebeforecombining' num2str(average_before_combining_terms) ...
            '-varcenter' num2str(variance_centering) ...
            '-' correction_method '-nfolds' num2str(n_folds)...
            '-' regularization_metric];
        box off;
        
        if noise_corrected_crossval
            fname = [fname '-noisecorrCV'];
        end

        fname_full_path = [figure_directory '/' fname];
        set(gcf, 'PaperSize', [8 8]);
        set(gcf, 'PaperPosition', [0.25 0.25 7.5 7.5]);
        print([fname_full_path '.pdf'], '-dpdf');
        % export_fig([figure_directory '/' fname '.pdf'], '-transparent', '-pdf');
        % export_fig([figure_directory '/' fname '.png'], '-png', '-r150');
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
