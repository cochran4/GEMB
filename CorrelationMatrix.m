function R = CorrelationMatrix(name)


% Load file
X  = readtable(name,'FileType','text');
ng = height(X)-1;

% Build table
T = table('Size',[ng,10],'VariableTypes',...
           {'string','string','double','double',...
           'double','double','double','double',...
           'double','cell'},...
'VariableNames',{'GENE','CHR','START','STOP','NSNPS',....
                           'NPARAM','N','VAL','ZSTAT','R'});
                     
% Loop through genes
for g=2:height(X)
    
    newStr        = split(X{g,:}{1});
    T.GENE(g-1)   = newStr{1};
    T.CHR(g-1)    = newStr{2};
    T.START(g-1)  = str2double(newStr{3}); 
    T.STOP(g-1)   = str2double(newStr{4}); 
    T.NSNPS(g-1)  = str2double(newStr{5});
    T.NPARAM(g-1) = str2double(newStr{6});
    T.N(g-1)      = str2double(newStr{7});
    T.VAL(g-1)    = str2double(newStr{8});
    T.ZSTAT(g-1)  = str2double(newStr{9});
    
    tmp = [];
    for j=10:length(newStr)
        tmp(j-9) = str2double(newStr{j}); 
    end
    T.R(g-1) = { tmp };
    
end

% Adjust correlation matrix to be positive definite
% Number of genes
ng = height(T);

% Build correlation matrix
R    = speye(ng,ng);
adjR = sparse(ng,ng);

% Adjust for non-positive definiteness
chr = unique(T.CHR);
tol = 10^(-6);
for i=1:length(chr)

    % Perform cholesky factorization
    ix    = find( ismember(T.CHR,chr(i)) )';
    
    % Collect
    for g=ix
         rho = T.R{g};
         if ~isempty(rho)
             R(g,g-(length(rho):-1:1)  ) = rho;
             R(g-(length(rho):-1:1),  g) = rho;
         end
    end
    
    % Cholesky factorization
    [~,p] = chol( R(ix,ix) );

    % Adjust negative eigenvalues
    if p~=0
        [V,D]    = eig(full(R(ix,ix)));
        D        = max(diag(D),tol);
        R(ix,ix) = (V*diag(D))/V; 
    end

end

