function IntracellularCalcium(type)
% Run gene-set test using weights from Ashhad-Narayanan model

% Add GEMB to path
addpath('..\.')

% Get parameter weights
load('AshhadNarayanan2013\Weights')

% Get calcium information 
Ca = readtable('CalciumGenes.xlsx');

% Get gene set names
GeneSetNames = Ca.Gene;

% Initialize weights
GeneSetWeights = zeros(height(Ca),1);

% Loop through parameters
for i=1:length(weights)
    GeneSetWeights( Ca.Index == i ) = weights(i)/sum( Ca.Index == i );
end

% Remove weights that are zero
GeneSetNames   = GeneSetNames(   GeneSetWeights~=0 );
GeneSetWeights = GeneSetWeights( GeneSetWeights~=0 );

% Normalize weights
GeneSetWeights = abs(GeneSetWeights)/sum(GeneSetWeights);

% Run gene set test
[~,p1] = AnalyzeData(GeneSetWeights,GeneSetNames,type,'NCBI37.3.gene.loc.txt','magma.genes.out.txt');

%-----------------------------------------------------------------------
% Remove CACNAC1


% Remove CACNA1C
Ca( ismember(Ca.Gene,'CACNA1C'), :) = [];

% Get gene set names
GeneSetNames = Ca.Gene;

% Initialize weights
GeneSetWeights = zeros(height(Ca),1);


% Loop through parameters
for i=1:length(weights)
    GeneSetWeights( Ca.Index == i ) = weights(i)/sum( Ca.Index == i );
end

% Remove weights that are zero
GeneSetNames   = GeneSetNames(   GeneSetWeights~=0 );
GeneSetWeights = GeneSetWeights( GeneSetWeights~=0 );

% Normalize weights
GeneSetWeights = abs(GeneSetWeights)/sum(GeneSetWeights);

% Run analysis
[~,p6] = AnalyzeData(GeneSetWeights,GeneSetNames,type,'NCBI37.3.gene.loc.txt','magma.genes.out.txt');

%------------------------------------
% Compare with KEGG Calcium signaling
load('CalciumGenes')    
 
% Assume uniform weights
CalciumWeights = ones(length(CalciumGenes),1);
CalciumWeights = CalciumWeights/sum(CalciumWeights);
 
[~,p2] = AnalyzeData(CalciumWeights,CalciumGenes,type,'NCBI37.3.gene.loc.txt','magma.genes.out.txt');
 
% Apply over-representation analysis
p3     = ORA(CalciumGenes,'NCBI37.3.gene.loc.txt','magma.genes.out.txt');

%------------------------------------
% Compare with CACNAC
load('CalciumGenes')    
[~,p4] = AnalyzeData(1,{'CACNA1C'},'Monte Carlo','NCBI37.3.gene.loc.txt','magma.genes.out.txt');
p5     = ORA({'CACNA1C'},'NCBI37.3.gene.loc.txt','magma.genes.out.txt');

% Output results
disp(['Ashhad-Narayanan:             ',sprintf('%0.3e',p1)])
disp(['KEGG Calcium:                 ',sprintf('%0.3f',p2)])
disp(['KEGG Calcium (ORA):           ',sprintf('%0.3f',p3)])
disp(['CACNA1C:                      ',sprintf('%0.3e',p4)])
disp(['CACNA1C (ORA):                ',sprintf('%0.3f',p5)])
disp(['Ashhad-Narayanan w/o Calcium: ',sprintf('%0.3e',p6)])
