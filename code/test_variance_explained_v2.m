
alpha = linspace(1,3,10);

z = 3;
T = 100;
b = randn(T,1);
fnonlin = @(f)(sin(f*alpha(z)) * b);

%% Estimate the correlation between the true nonlinear function and the best linear approximation

methods = {'least-squares','ridge'};
K = {[],2.^(-60:60)};

true_r = nan(length(N),length(methods));
for k = 1:length(methods);
    N = round(linspace(30,300,5));
    for i = 1:length(N)
        fprintf('%s,%d\n', methods{k}, N(i));
        drawnow;
        n_smps = 10;
        r = nan(n_smps,1);
        for j = 1:n_smps
            % true correlation
            F = randn(N(i),T);
            y = fnonlin(F);
            prediction = regress_predictions_from_3way_crossval(F, y, 10, methods{k}, K{k});
            r(j) = corr(prediction,y);
        end
        true_r(i,k) = mean(r);
    end
end


%% Estimate the correlation between the true nonlinear function and the best linear approximation

N = 1000;
n_smps = 10;
r = nan(n_smps,1);
for j = 1:n_smps
    % true correlation
    F = randn(N,T);
    y = fnonlin(F);
    r(j) = corr(fnonlin(F), F*pinv(F)*fnonlin(F));
    
end
true_r = mean(r)



%%

root_directory = '/Users/svnh2/Desktop/projects';
addpath([root_directory '/general-analysis-code']);

N = 100;
n_smps = 10;
r = nan(n_smps,1);
for i = 1:n_smps
    F = randn(N,T);
    y1 = fnonlin(F) + randn(N,1);
    y2 = fnonlin(F) + randn(N,1);
    prediction1 = regress_predictions_from_3way_crossval(F, y1, 10, 'ridge');
    prediction2 = regress_predictions_from_3way_crossval(F, y2, 10, 'ridge');
    
    r1 = corr(y1, prediction1);
    r2 = corr(y2, prediction2);
    if corr(y1,y2) > 0 && corr(prediction1,prediction2) > 0
        r(i) = (r1/2 + r2/2) / sqrt(corr(y1,y2) * corr(prediction1,prediction2));
    end
end

estimated_r = nanmedian(r)
