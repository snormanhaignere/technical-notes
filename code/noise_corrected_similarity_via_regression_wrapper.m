correction_methods = {'variance-based', 'correlation-based'};
average_before_combining_terms = [true, false];

for i = 1:length(correction_methods)
    for j = 1:length(average_before_combining_terms)
        switch correction_methods{i}
            case 'variance-based'
                variance_centering = [true false];
            case 'correlation-based'
                variance_centering = false;
            otherwise
                error('Switch statement fell through');
        end
        
        for k = 1:length(variance_centering)
            noise_corrected_similarity_via_regression(...
                correction_methods{i}, average_before_combining_terms(j), ...
                variance_centering(k))
        end
    end
end