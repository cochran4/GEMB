function RunPRCC
% Save partial rank correlation coefficients for each parameter

% Load parameter design matrix
load('X')

% Collect data in a vector
for k=1:size(X,1)
    tmp          = load(['ave_cai',num2str(k),'.dat']);
    ave_cai(k,1) = tmp(1);
end

% Calculate PRCC
rho     = partialcorr( [X(:,iy),ave_cai(:,1)],'Type','Spearman');
weights = rho(1:end-1,end);

% Save weights
save('Weights','weights')