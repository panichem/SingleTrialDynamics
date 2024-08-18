# SingleTrialDynamics

This repository contains the custom behavioral and analysis functions supporting Panichello, Jonikaitis, Oh, Zhu, Trepka, & Moore, 2024.

./task/taskCode.m presents stimuli and record behavioral responses in the memory-guided saccade task. 

./clust/clusterMassOneSampZ.m is the core function used to identify significant epochs of above-chance confidence ('On states'), corrected for multiple comparisons, as described in Figure 4. 

./bmm/bmmFit.m and ./bmm/bmmFit_1 return the log likehood and parameters of the best-fitting double- and single-component beta mixture models for data with domain [0, 1], as described in Figure ED5. 

./ccgs/analyze_ccgs.m is the core function used to compute the firing-rate normalized and jitter corrected CCGs analyzed in Figures 5 and 6. 
