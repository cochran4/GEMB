# GEMB
Gene-set enrichment with math biology

This code accompanies the paper "Gene-Set Enrichment with Mathematical Biology" by AL Cochran, DB Forger, S Zoellner, and MG McInnis. Please cite the paper appropriately. 

The code requires Matlab to run.  The main function is "WeightedGSTest" which runs a weighted gene-set test to recover a p-value as discussed in the paper above.  For example, suppose that you calculated a weighted gene-set statistic of v=1000, there were n=10,000 genes analyzed, and gene weights were [0.1,0.1,0.3,0.4,0.1]. Then, you can estimate a p-value for the statistic by running the following command in matlab:

  p = WeightedGSTest(1000,[0.1,0.1,0.3,0.4,0.1],10000,'Monte Carlo');

using a Monte Carlo approximation. Alternatively, a quicker estimate could be achieved using 

  p = WeightedGSTest(1000,[0.1,0.1,0.3,0.4,0.1],10000,'Normal');

The other functions require that you gene weights, gene names, and two files containing gene associations and gene locations. The latter files are not available in this repository, because they use data from other sources (e.g., data from the Psychiatry Genetics Consortium, PCG). We used PCG data,  


