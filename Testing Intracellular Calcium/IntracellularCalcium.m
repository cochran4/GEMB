function IntracellularCalcium(method,rankfile,locfile,corrfile)
% Run gene-set test using weights from Ashhad-Narayanan model
% method: 'Monte Carlo' or 'Normal' or 'Correlation'
% rankfile: file name where gene ranks are stored
% locfile:  file name where gene locations are stored
% corrfile: file name where gene correlations are stored

% Add GEMB to path
addpath('..\')

% Get parameter weights
Ca = readtable('Weights_Calcium.xlsx');

% Get gene set names
GeneSetNames = Ca.Gene;

% Initialize weights
GeneSetWeights = Ca.Weight;

% Run gene set test
if nargin > 3
    [~,p1] = AnalyzeData(GeneSetWeights,GeneSetNames,method,locfile,rankfile,corrfile);
else
    [~,p1] = AnalyzeData(GeneSetWeights,GeneSetNames,method,locfile,rankfile);
end
%-----------------------------------------------------------------------
% Remove CACNAC1

% Remove CACNA1C
Ca( ismember(Ca.Gene,'CACNA1C'), :) = [];

% Get gene set names
GeneSetNames = Ca.Gene;

% Initialize weights
GeneSetWeights = Ca.Weight;

% Normalize weights
GeneSetWeights = abs(GeneSetWeights)/sum(GeneSetWeights);

% Run analysis
if nargin > 3
    [~,p6] = AnalyzeData(GeneSetWeights,GeneSetNames,method,locfile,rankfile,corrfile);
else
    [~,p6] = AnalyzeData(GeneSetWeights,GeneSetNames,method,locfile,rankfile);
end

%------------------------------------
% Compare with KEGG Calcium signaling

% Get genes in KEGG calcium signaling
T = readtable('KEGG_Calcium_Genes.xlsx');

% Names
KEGGCalciumGenes = T.Name;
 
% Assume uniform weights
KEGGCalciumWeights = ones(length(KEGGCalciumGenes),1);
KEGGCalciumWeights = KEGGCalciumWeights/sum(KEGGCalciumWeights);
 
if nargin > 3
    [~,p2] = AnalyzeData(KEGGCalciumWeights,KEGGCalciumGenes,method,locfile,rankfile,corrfile); 
else
    [~,p2] = AnalyzeData(KEGGCalciumWeights,KEGGCalciumGenes,method,locfile,rankfile);
end
 

% Apply over-representation analysis
p3     = ORA(KEGGCalciumGenes,locfile,rankfile);

%------------------------------------
% Compare with testing with only CACNAC

[~,p4] = AnalyzeData(1,{'CACNA1C'},method,locfile,rankfile);
p5     = ORA({'CACNA1C'},locfile,rankfile);

% Output results
disp(['Ashhad-Narayanan:             ',sprintf('%0.3e',p1)])
disp(['KEGG Calcium:                 ',sprintf('%0.3f',p2)])
disp(['KEGG Calcium (ORA):           ',sprintf('%0.3f',p3)])
disp(['KEGG CACNA1C:                 ',sprintf('%0.3e',p4)])
disp(['KEGG CACNA1C (ORA):           ',sprintf('%0.3f',p5)])
disp(['Ashhad-Narayanan w/o Calcium: ',sprintf('%0.3e',p6)])
