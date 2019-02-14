function p = ORA(GeneSetNames,LocFile,RankFile)
% Over-representation analysis
% Input:  
%           GeneSetNames   = m-by-1 cell of names of gene in gene set
%           LocFile        = file name of gene location data
%           RankFile       = file name of gene rank data

% Load gene location data
GeneLoc   = readtable(LocFile);

% Load gene rank data
GeneRank  = readtable(RankFile);

%--------------------------------------------------------------------------
% Find id of gene names
%--------------------------------------------------------------------------

% Find gene in gene location file
[~,ix] = ismember(GeneSetNames,GeneLoc.Var6);

% Warning: some names in gene set are unavailable
if any(ix==0)
   disp(['Warning: the following gene names are not in the gene location file....'])
   missing = find(ix==0);
   for i=1:length(missing)
        disp(['     ',GeneSetNames{missing(i)}])
   end
   disp('Removing these genes from subsequent analysis.')
end

% Remove genes that are not available
GeneSetNames   = GeneSetNames(   ix~=0 );
ix             = ix( ix~=0 );

% Find gene set indices
GeneSetID      = GeneLoc.Var1(ix);
[~,GeneSetIX]  = ismember( GeneSetID, GeneRank.GENE ); 

% Warning: some names in gene set are unavailable
if any(GeneSetIX==0)
   disp(['Warning: the following gene names are not in the gene rank file....'])
   missing = find(GeneSetIX==0);
   for i=1:length(missing)
        disp(['     ',GeneSetNames{missing(i)}])
   end
   disp('Removing these genes from subsequent analysis.')
end

% Remove genes that are not available
GeneSetNames   = GeneSetNames(   GeneSetIX~=0 );
GeneSetIX      = GeneSetIX( GeneSetIX~=0 );

% Sort P-values from gene-level analysis
P = GeneRank.P;
disp(['Number of genes: ',num2str(length(P))])
[Psorted,IX]   = sort( P, 'ascend' ); 
n        = length(P);

% False discovery rate FDR
k = find( Psorted(:) <= 0.1*(1:n)' / n , 1, 'last' );

% Adjust if not enough significant genes based on FDR
if isempty(k) || k < 200
    k = round( 0.01*n ); % top 1% of genes
end

% Number of genes in gene set
m = length(GeneSetIX);

% Check for genes in gene-set
[~,GeneSetRanks] = ismember( GeneSetIX , IX ); 

% Over-represented genes in gene set
or = sum( GeneSetRanks  <= k );

% Perform one-sided hypergeometric
p = 1-hygecdf(or-1,n,k,m); 