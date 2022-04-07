# Naive-VS-Efficient-Clustering-Comparison-LFR
Code for comparing naive and efficient undirected likelihood, modularity and Louvain Method clustering. Automatically generates, and tests objective functions on, LFR benchmark networks.

This code was created to support the paper "Improving computational performance in likelihood-based network clustering" by Alan Ballard and Marcus B. Perry (Stat. 2019; 8:e256. https://doi.org/10.1002/sta4.256). 

This code will generate a set of the so-called LFR benchmark undirected networks proposed by Lancichinetti et al (2008) in "Benchmark graphs for testing community detection algorithms", along with corresponding ground-truth solutions.
Each LFR network is then clustered into a range of user-specified number of clusters using the likelihood objective function proposed by Perry et al (2013) in "On the statistical detection of clusters in undirected networks", the efficient version of this likelihood objective function proposed in the current paper and the modularity objective function proposed by Newman (2006) in "Modularity and community structure in networks". 
This clustering is accomplished using simulated annealing and the LFR network generation parameters, the cluster number range, and the simulated annealing cooling schedule are adjustable within the code.

Each network is also clustered using the so-called Louvain Method of optimizing modularity proposed by Blondel et al (2008) in "Fast unfolding of communities in large networks". This is an agglommerative, greedy clustering method that is known for its speed and although it differs from simulated annealing in optimization mechanics, it is included here for comparison of clustering quality and time-to-solution.
