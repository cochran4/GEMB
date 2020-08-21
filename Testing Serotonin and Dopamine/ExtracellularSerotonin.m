function ExtracellularSerotonin(method,rankfile,locfile,corrfile)
% Run gene-set test using weights from Ashhad-Narayanan model
% method: 'Monte Carlo' or 'Normal' or 'Correlation'
% rankfile: file name where gene ranks are stored
% locfile:  file name where gene locations are stored
% corrfile: file name where gene correlations are stored

% Add GEMB to path
addpath('..\')

% Get parameter weights
Da = readtable('Weights_Serotonin.xlsx');

% Get gene set names
GeneSetNames = Da.Gene;

% Initialize weights
GeneSetWeights = Da.Weight;

% Run gene set test
if nargin > 3
    [~,p1] = AnalyzeData(GeneSetWeights,GeneSetNames,method,locfile,rankfile,corrfile);
else
    [~,p1] = AnalyzeData(GeneSetWeights,GeneSetNames,method,locfile,rankfile);
end


% Output results
disp(['Serotonin model:             ',sprintf('%0.3f',p1)])
