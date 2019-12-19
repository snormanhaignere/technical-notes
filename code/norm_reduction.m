% Two equations that show how to compute the reliability metric used in the
% Neuron paper. The formula in the paper is inaccurate because all of the L2
% norms should be squared.

v1 = rand(10,1);
v2 = rand(10,1);

e = v1 - v2 * pinv(v2) * v1;
r = 1 - sum(e.^2)/sum(v1.^2)

proj_v1 = v2 * (v2'/(norm(v2).^2)) * v1;
e = v1 - proj_v1;
r = 1 - (norm(e).^2)/(norm(v1).^2)