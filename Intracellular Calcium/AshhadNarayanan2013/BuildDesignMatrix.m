
% Latin hypercube sampling
rng(0)
X = 1 + lhsnorm(zeros(13,1),eye(13)*0.05^2,320);
save('X','X')
