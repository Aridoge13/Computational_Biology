

Using BFGS method to optimize GTR rate parameters, to disable this specify "--no-bfgs" 



This is RAxML version 8.2.12 released by Alexandros Stamatakis on May 2018.

With greatly appreciated code contributions by:
Andre Aberer      (HITS)
Simon Berger      (HITS)
Alexey Kozlov     (HITS)
Kassian Kobert    (HITS)
David Dao         (KIT and HITS)
Sarah Lutteropp   (KIT and HITS)
Nick Pattengale   (Sandia)
Wayne Pfeiffer    (SDSC)
Akifumi S. Tanabe (NRIFS)
Charlie Taylor    (UF)


Alignment has 564 distinct alignment patterns

Proportion of gaps and completely undetermined characters in this alignment: 10.47%

RAxML rapid hill-climbing mode

Using 1 distinct models/data partitions with joint branch length optimization


Executing 1 inferences on the original alignment using 1 distinct randomized MP trees

All free model parameters will be estimated by RAxML
GAMMA model of rate heterogeneity, ML estimate of alpha-parameter

GAMMA Model parameters will be estimated up to an accuracy of 0.1000000000 Log Likelihood units

Partition: 0
Alignment Patterns: 564
Name: No Name Provided
DataType: DNA
Substitution Matrix: GTR




RAxML was called as follows:

/usr/bin/raxmlHPC-PTHREADS-AVX -T 4 -s /home/jijo/Projects/Comp_Bio/resources/alignment/alignment.fasta -n bison -m GTRGAMMA -p 12345 -w /home/jijo/Projects/Comp_Bio/resources/tree 


Partition: 0 with name: No Name Provided
Base frequencies: 0.307 0.266 0.125 0.302 

Inference[0]: Time 2.079827 GAMMA-based likelihood -7318.050632, best rearrangement setting 5
alpha[0]: 0.489262 rates[0] ac ag at cg ct gt: 5.130684 7.406115 2.048941 0.532938 14.833458 1.000000 


Conducting final model optimizations on all 1 trees under GAMMA-based models ....

Inference[0] final GAMMA-based Likelihood: -7317.892769 tree written to file /home/jijo/Projects/Comp_Bio/resources/tree/RAxML_result.bison


Starting final GAMMA-based thorough Optimization on tree 0 likelihood -7317.892769 .... 

Final GAMMA-based Score of best tree -7317.892769

Program execution info written to /home/jijo/Projects/Comp_Bio/resources/tree/RAxML_info.bison
Best-scoring ML tree written to: /home/jijo/Projects/Comp_Bio/resources/tree/RAxML_bestTree.bison

Overall execution time: 2.488230 secs or 0.000691 hours or 0.000029 days

