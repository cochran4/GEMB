function p = WeightedGSTest(v,w,n,approach)
% Calculate P-value for weighted gene-set test
% Input:
%       v        - weighted average rank 
%       w        - gene weights held in vector of length m
%       n        - number of genes
%       approach - 'Monte Carlo' or 'Normal'

% Intialize variables
par.m     = length(w);
par.n     = n;
w         = w(:); 
par.v     = v;
par.tol   = 10^(-6);                % tolerance

% Two approaches to estimate one-sided p-value
switch approach
    case 'Monte Carlo'
        par.minit = 10^7;                   % minimum iteration
        conv      = MonteCarloApproach(w,par); 
        p         = conv(end);
    case 'Normal'
        mu   = (n+1)/2;
        sig2 = sum(w(:).^2)*(n^2-1)/12;
        p    = normcdf( (v-mu)/sqrt(sig2) );
end

function p = MonteCarloApproach(w,par)
% Run monte-carlo until convergence

% Initialize variables
it   = 1;
cnt  = 0;
conv = zeros(par.minit,1);

% Repeatedly draw samples
while it <= par.minit || abs(conv(it-1)-conv(it-2))/conv(it-2) > par.tol

    % Draw random subset
    ix = randperm(par.n,par.m);

    % Estimate weighted average; count if more extreme than x            
    cnt  = cnt + (ix*w <= par.v);           

    % Running convergence
    conv(it,1) = cnt/it;

    % Update iteration
    it = it + 1; 

end    

% Get last estimate
p = conv(end);