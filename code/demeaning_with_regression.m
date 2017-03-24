X = rand(10,3);
y = rand(10,1);
b = pinv([ones(10,1), X]) * y
b = pinv(bsxfun(@minus, X, mean(X))) * (y - mean(y))