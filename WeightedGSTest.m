function [p,conv] = WeightedGSTest(v,w,n,approach,R,GeneSetIX)
% Calculate P-value for weighted gene-set test
% Input:
%       x         - weighted average rank 
%       w         - gene weights held in vector of length m
%       n         - number of genes
%       approach  - 'Monte Carlo' or 'Normal' or 'Correlation'
%       R         - gene-by-gene correlation matrix (required for 'Correlation')
%       GeneSetIX - gene set indices (required for 'Correlation')

% Intialize variables
par.m     = length(w);
par.n     = n;
w         = w(:); 
par.v     = v;
par.tol   = 10^(-6);                % tolerance

% Correlation matrix
if nargin > 3 && ~strcmp(approach,'Correlation')
    disp('Gene correlation matrix is only used for Monte Carlo approach. \n')
    disp('Ignoring gene correlation matrix. \n')
end

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
    case 'Correlation'
        par.tol   = 10^(-2);                % higher tolerance b/c of slower computation
        par.minit = 5*10^6;                 % lower min iteration b/c of slower computation
        conv      = CorrelationApproach(w,par,R,GeneSetIX); 
        p         = conv(end);
end

function conv = MonteCarloApproach(w,par)
% Run monte-carlo until convergence

% Initialize variables
it    = 0;
cnt   = 0;

while it <= par.minit || abs(conv(it-1)-conv(it-2))/conv(it-2) > par.tol

    % Draw random ranks
    ix = randperm(par.n,par.m);
    
    % Estimate weighted average; count if more extreme than x            
    cnt  = cnt + sum(ix'*w <= par.v);           

    % Update iteration
    it = it + batch;    
   
    % Running convergence
    conv(it,1) = cnt/it;

end       

function conv = CorrelationApproach(w,par,R,GeneSetIX)
% Run monte-carlo until convergence

% Initialize variables
it    = 0;
cnt   = 0;
V     = chol(R)';
batch = 10^4;


while it <= par.minit || abs(conv(it/batch)-conv(it/batch-1)) > par.tol

    % Draw random ranks
    parfor j=1:10
        iy{j} = tiedrank( V*randn(par.n,batch/10) );
    end
    
    ix = [];
    for j=1:10
        ix = [ix,iy{j}(GeneSetIX,:)];
    end
    
    % Estimate weighted average; count if more extreme than x            
    cnt  = cnt + sum(ix'*w <= par.v);           

    % Update iteration
    it = it + batch;    
   
    % Running convergence
    conv(it/batch,1) = cnt/it;
    
end       