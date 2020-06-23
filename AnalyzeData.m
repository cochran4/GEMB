function [v,p] = AnalyzeData(GeneSetWeights,GeneSetNames,Method,LocFile,RankFile,CorrFile)
% Apply weighted gene set analysis to given geneset
% Input:  
%           GeneSetNames   = m-by-1 cell of names of gene in gene set
%           GeneSetWeights = m-by-1 vector of respective weights in gene set 
%           Method         = 'Monte Carlo' or 'Normal'
%           LocFile        = file name of gene location data
%           RankFile       = file name of gene rank data
%           CorrFile       = file name of gene correlation data

% Load gene location data
GeneLoc   = readtable(LocFile,'FileType','text');

% Load gene rank data
GeneRank  = readtable(RankFile,'FileType','text');

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
GeneSetWeights = GeneSetWeights( ix~=0 ); 
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
GeneSetWeights = GeneSetWeights( GeneSetIX~=0 ); 
GeneSetIX      = GeneSetIX( GeneSetIX~=0 );

%--------------------------------------------------------------------------
% Check weights are non-negative and sum to one
%--------------------------------------------------------------------------
if any(GeneSetWeights<0)
    disp('Taking absolute value of weights...')
    GeneSetWeights = abs(GeneSetWeights);
end

if sum(GeneSetWeights) ~=0
    disp('Normalizing weights ...')
    GeneSetWeights = GeneSetWeights / sum( GeneSetWeights );
end


%--------------------------------------------------------------------------

% Sort P-values from gene-level analysis
P = GeneRank.P;
disp(['Number of genes: ',num2str(length(P))])
disp([])
[~,IX]   = sort( P, 'ascend' ); 
n        = length(P);

% Get ranks of genes in gene set
[~,GeneSetRanks] = ismember( GeneSetIX , IX ); 
v                = GeneSetRanks(:)'*GeneSetWeights(:);  % Statistic
    

% Add in correlation matrix
if nargin == 5
    p                = WeightedGSTest(v,GeneSetWeights,n,Method); 
else
    % Load correlation data
    R = CorrelationMatrix(CorrFile);
    p = WeightedGSTest(v,GeneSetWeights,n,Method,R,GeneSetIX); 
end
