# GEMB
Gene-set enrichment with math biology

This code accompanies the paper "Gene-Set Enrichment with Mathematical Biology" by AL Cochran, K Nieser, DB Forger, S Zoellner, and MG McInnis at GigaScience. Please cite the paper appropriately. 

The code requires Matlab to run.  The main function is "WeightedGSTest" which runs a weighted gene-set test to recover a p-value as discussed in the paper above.  For example, suppose that you calculated a weighted gene-set statistic of v=1000, there were n=10,000 genes analyzed, and weights for genes in your target set were [0.1,0.1,0.3,0.4,0.1]. Then, you can estimate a p-value for the statistic by running the following command in matlab:

  p = WeightedGSTest(1000,[0.1,0.1,0.3,0.4,0.1],10000,'Monte Carlo');

using a Monte Carlo approximation. Alternatively, a quicker but less accurate estimate could be achieved using 

  p = WeightedGSTest(1000,[0.1,0.1,0.3,0.4,0.1],10000,'Normal');

The other functions require that you gene weights, gene names, and two files containing gene associations and gene locations. The latter files are not available in this repository, because they use data from other sources (e.g., data from the Psychiatry Genetics Consortium, PCG). We were able to construct these files with publically-available data and gene-set analysis software magma. If you were able to recover such files, which we labeled 'NCBI37.3.gene.loc.txt' (gene locations) and 'magma.genes.out.txt' (gene associations), you could run the command in matlab to run the gene-set test using weights in GeneSetWeights

  [v,p] = AnalyzeData(GeneSetWeights,GeneSetNames,type,'NCBI37.3.gene.loc.txt','magma.genes.out.txt');
  
or run a typical overrepresentation analysis with

  p     = ORA(CalciumGenes,'NCBI37.3.gene.loc.txt','magma.genes.out.txt');


The calcium example in the association paper can be reproduced with files such as 'NCBI37.3.gene.loc.txt' (gene locations) and 'magma.genes.out.txt' (gene associations) by running

  IntracellularCalcium('Monte Carlo')
  
or

  IntracellularCalcium('Normal')
  
depending on what method you want to use to estimate p values. 

Finally, you can find associated genes and PRCC values based on intracellular calcium ion concentrations predicted by the Ashhad & Narayanan model in 'IntracellularCalcium\CalciumGenes.xlsx' and 'IntracellularCalcium\AshhadNarayanan2013\Weights.mat'. Results and certain neuron code for running the Ashhad & Naraynan model are in the folder 'IntracellularCalcium\AshhadNarayanan2013', but all files from the original Ashhad & Naraynan model found on modelDB have been removed. 



