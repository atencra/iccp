Inferior Colliculus Neuron Pairs Analysis
===========

These functions analyze and plot data from paired neuron analysis. The data was recorded from the auditory midbrain. The specific site was the central nucleus of the Inferior Colliculus. 

The recordings were made with multi-channel electrodes, where each channel was separated by 150 microns, and the channels were linearly organized.

The goal of the project was to compare the processing of neighboring neurons in the auditory neurons. These neurons were recorded from the same electrode channel, and thus they are located within 75 um of each other. To extract the neurons, the recording trace from each channel must first be spike-sorted. Next, I estimated the cross-correlation from each pair of neuorns from an electrode contact. Positive peaks around zero delay in the function indicate synchronized responses and functional connectivity. Next, how the pairs of neurons process sounds was compared. Finally, I also estimated information encoding using information theory techniques.

Each file begins with iccp, since that is shorthand for inferior colliculus pairs. 

Figures and statistical analyses using these functions have been included in a paper that was published at the journal Neuroscience.

You can read the paper at: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5031551/

