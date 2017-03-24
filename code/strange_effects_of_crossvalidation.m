% construct basis functions with decaying eigenvalues
N = 100;
T = 10:5:200;
n_T = length(T);

n_iter = 100;
rU = nan(n_iter, n_T);
rF = nan(n_iter, n_T);
for j = 1:n_T
    for i = 1:100
        %         [U,~,V] = svd(randn(N,T(j)),'econ');
        %         S = diag(1./(1:min(N,T(j))));
        %         F = U*S*V';
        %
        %         % create a signal that is a mixture of the first P principal components
        %         P = T(j);
        %         sig = U(:,1:P) * randn(P,1);
        %
        
        F = randn(N,T(j));
        [U,S,V] = svd(F,'econ');
        %         S = eye(min(N,T(j)));
        %         S = diag(1./(1:min(N,T(j))));
        % S = diag((1:min(N,T(j))));
        F = U*S*V';
        sig = F * rand(T(j),1);
        
        % add noise
        noise = 0*sqrt(var(sig)/2) * randn(N,1);
        data = sig + noise;
        
        % first half for training
        % second half for testing
        train_inds = 1:N/2;
        test_inds = N/2+1:N;
        
        % predict responses using U
        % predictions are poor
        bh = pinv(U(train_inds,:)) * data(train_inds);
        rU(i,j) = corr(data(test_inds), U(test_inds,:) * bh);
        
        % predict responses using X
        % predictions are much better, why?
        bh = pinv(F(train_inds,:)) * data(train_inds);
        rF(i,j) = corr(data(test_inds), F(test_inds,:) * bh);
        
    end
end
figure;
plot(T, [median(rF)', median(rU)'], 'LineWidth', 2);
legend('X', 'U');
xlim([0 max(T)]);
yL = ylim;
hold on;
plot([N/2, N/2], ylim, 'k--', 'LineWidth', 2);
legend('F', 'U', 'N/2');
xlabel('Number of Features');
ylabel('Cross-validated Correlation');
set(gca, 'FontSize', 16);

%%
% set(gcf, 'PaperSize', [6 6]);
% set(gcf, 'PaperPosition', [0.25 0.25 5.5 5.5]);
% print([pwd '/strange_effects_of_crossval.pdf'],'-dpdf')
% print([pwd '/strange_effects_of_crossval.png'],'-dpng', '-r200');
%
