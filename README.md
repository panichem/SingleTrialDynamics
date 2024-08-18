# SingleTrialDynamics

This repository contains core behavioral and analysis functions supporting Panichello, Jonikaitis, Oh, Zhu, Trepka, & Moore, 2024.

taskCode.m is the script used to present stimuli and record behavioral responses in the memory-guided saccade task. 

logregTFdir is the core function used to train and test the classifiers presented in Figure 2.

computePhaseModulationIndex.m and the helper function phaseModulationIndex.m are the core functions used to compute the z-scored color information statistic presented in Figure 3. clusterMassOneSamp.m, clusterMassDependent.m, and the helper function getClust.m are the core functions used compute the significance of this color information statistic using cluster-corrected t-tests.

planesDemo.m with the helper function planeAngle.m loads example firing rate data from planesDemo.mat, visualizes population firing rates for each condition in a low-D space, and computes the angle between the upper and lower color planes, as described in Figure 4.
