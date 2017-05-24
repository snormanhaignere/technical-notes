correction_methods = {'variance-based'};
average_before_combining_terms = true;
n_folds = [2, 4];
n_features = 10;
noise_factor = [2, 3];
noise_corrected_crossval = [true, false];
regularization_metric = {'demeaned-squared-error', 'unnormalized-squared-error', 'pearson'};

for n = 1:length(noise_factor)
    for m = 1:length(n_features)
        for l = 1:length(n_folds)
            for i = 1:length(correction_methods)
                for j = 1:length(average_before_combining_terms)
                    switch correction_methods{i}
                        case 'variance-based'
                            variance_centering = false;
                        case 'correlation-based'
                            variance_centering = false;
                        otherwise
                            error('Switch statement fell through');
                    end
                    
                    for k = 1:length(variance_centering)
                        for o = 1:length(noise_corrected_crossval)
                            for p = 1:length(regularization_metric)
                                noise_corrected_similarity_via_regression(...
                                    correction_methods{i}, average_before_combining_terms(j), ...
                                    variance_centering(k), n_folds(l), n_features(m), ...
                                    noise_factor(n), noise_corrected_crossval(o), regularization_metric{p});
                            end
                        end
                    end
                end
            end
        end
    end
end

