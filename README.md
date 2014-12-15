Multi-Edge Randomizer: Maximum entropy Multi-Edge network generator
========================================================================

 Copyright 2013 Oleguer Sagarra and Francesc Font-Clos. All rights reserved. Code under License GPLv3.
______________________________________________________________________________________

The present program may be used to generate maximally random multi-edge networks with prescribed strength sequences using 3 different ensembles: Micro-Canonical (rewirings), Canonical or Grand-Canonical.
Multi-edge networks are weighted networks composed by distinguishable weights (weights which may be univocally identified such as time-stamped interactions or users riding on transportation systems).

Randomizing a network keeping some properties of the original network can be useful for studying the effect of a certain feature on any topological or dynamical property of such network.
Using different ensembles may be useful for performance issues: For the case of multi-edge networks, it can be proved analytically that all ensembles are strictly equivalent in the large event limit. 

Comparing real networks with their randomized counterparts is very useful for data analysis and detection of statistically significant patterns. With this software you migh do just that, but if you use our software please do cite us.

## References 

[1] Statistical mechanics of multiedge networks
	Sagarra O., Pérez-Vicente C. and Díaz-Guilera, A.  Phys. Rev. E 88, 062806 (2013)
    [Physical Review E](http://pre.aps.org/abstract/PRE/v88/i6/e062806)

[2] The configuration multi-edge model: Assessing the effect of fixing node strengths on weighted network magnitudes
	Sagarra O., Font-Clos F., Pérez-Vicente C. and Díaz-Guilera, A.
	[Europhysics Letters](http://iopscience.iop.org/0295-5075/107/3/38002/article;jsessionid=D22CAAF312F43653DA0C1279853CBF0C.c3)



## Requirements and Installation

  You only need to install the GSL library: http://www.gnu.org/software/gsl/
  
  How to tutorial: http://www.brianomeara.info/tutorials/brownie/gsl

## Compilation

  Simply run:

    $ make 

  The executable created is called **MultiEdgeGen**


## Execution

The executable is called **MultiEdgeGen**. The program has a set of options. To introduce on option you must introduce first its option key (preceded by a dash), one white space and then the argument value. The only compulsory option is the input network file (-f) in weighted adjacency list format, the number of nodes (-N) and the directed/undirected option (-d).
Arguments can appear in any order. If an argument does not appear the program gets the default value:

Examples
______________
 Generate ensemble averages for the maximum entropy grand-canonical network with prescribed strength sequence given by netFilename over 1000 reps and print an instance of such ensemble.
 	$ ./MutliEdgeRand -f netFILEname -d 1 -N 940 -r 1000 -p 1
 
 
## Options and Arguments

#### Options
```
	Compulsory:
		-N      	Number of nodes [int]
		-d			Directed (1) or undirected (0) option [int]
		-f file_s 	Path to file with strength sequence in form on each line: node_num(int) s_out(int) s_in(int) in the directed case, node_num (int) s(int) otherwise
 	Optional
		-s 			Initial seed for random generator [int] (default=1)
 		-e 			Ensemble_option: 0 (canonical, multinomial), 1 (grand-canonical, poisson + multinomial), 2 (grand-canonical, poisson indep.), 3(micro-canonical). [int] (default=2)
		-p 			Print adj list? 0 (no), 1 (yes) [int] (default=0) 
 		-x 			Exponent for log-binning on network stats (-1 for no log binning) [int] (default=-1)
 		-r 			Number of reps for averaging [int] (default=100)
 		-v 			Verbose (1 for on, 0 for off) [int] (default 0)
 		-c 			Clustering option (1 for yes) (warning: Depending on av_s makes simulations orders of magnitude slower) [int] (default=0)
 		-l 			Self-loop option (>0 for accepting them) (default =1)
        -w 			Compute analytic distribution of weights? (>0 for yes, takes some time) (default=0)
        -h 			Number of header lines on file_s [int] (default=1)
```

## Brief description of the Multi-Edge Generator

This program generates ensemble expectations of various properties for networks belonging to the ensemble of all maximum entropy multi-edge networks with prescribed strength sequence for 3 different types of ensembles.

In case you are interested in randomizing binary networks take a look at [RandNetGen](https://github.com/polcolomer/RandNetGen) by Pol Colomer-de-Simón.

## Outputs

For both directed and undirected cases:

The generated outputs are files with extension *.hist* if they correspond to network statistics (binned), *.tr* for a network instance of the ensemble in weighted adjacency format and *.list* for the node list with different features.


The files with the suffix *ens* reffer to averages over ensembles. The ones without the *ens* suffix refer to a single realization (and hence the *.list* files in such case do not contain averages nor std, but only raw values for each node).


**.list files**
   The file *ens_rXnode_list.list* contains ensemble averages for different node features, together with their standard deviations.
   For the undirected case: (Note that the clustering only appears if -c flag is set)
   
```
Node_num k k_std k_analitic k_analitic_std s s_std Y2 Y2_std k_nn k_nn_std k^w_nn k^w_nn_std s^w_nn s^w_nn_std (optionally  clust c_std  c^w c^w_std)
```
Explicit formulas: (in latex)
```
	k_i = \sum \Theta(t_ij)
	k^{anal}_i = \sum (1-e^{-s_is_j/T})
	s_i = \sum t_{ij}
	Y2 = \sum t_ij^2 / s_i
	k_nn_i = \sum \Theta(t_{ij}) k_j / k_i
	k^w_nn_i = \sum t_{ij} k_j /s_i
	s^w_nn = \sum t_{ij} s_j / s_i
	c = \sum_{jk} \Theta(t_ij)\Thetat(t_jk)\Theta(t_ki) / k_i(k_i-1)
	c^w = \sum_{jk} (t_ij + t_ik) \Theta(t_ij) \Thetat(t_jk)\Theta(t_ki) / [2s_i(k_i-1)]
```

For the directed case: (In this case the clustering is not defined)

```
Node_num k k_std k_analitic k_analitic_std s s_std Y2 Y2_std k_nn k_nn_std k^w_nn k^w_nn_std s^w_nn s^w_nn_std 
```
First for the out case, then for the in case. 

Explicit formulas (in latex):
Out:
```
	k_i^{out} = \sum_j \Theta(t_ij)
	k^{anal}_i^{out} = \sum_j (1-e^{-s_i^out}s_j^{in}/T})
	s_i^{out} = \sum_j t_{ij}, 
	Y2^{out} = \sum_j t_ij^2 / s_i^{out}
	k_nn_i^{out} = \sum_j \Theta(t_{ij}) k_j^{in} / k_i^{out}
	k^w_nn_i = \sum_j t_{ij} k_j^{in} /s_i^{out}, 
	s^w_nn = \sum_j t_{ij} s_j^{in} / s_i^{out}
In:
	k_j^{in} = \sum_j \Theta(t_ij)
	k^{anal}_j = \sum_j (1-e^{-s_i^out s_j^in/T})
	s_j^in = \sum_i t_{ij}
	Y2^{in} = \sum_i t_ij^2 / s_j^{in}
	k_nn^{in} = \sum_i \Theta(t_{ij}) k_i^{out} / k_j^{in}
	k^w_nn = \sum_i t_{ij} k_i^{out} /s_j^{in}, 
	s^w_nn^{in} = \sum_i t_{ij} s_i^{out} / s_j^{in}
```
Please note that the averages for Y2, Clustering and Average neighbor properties are performed only over realizations where the nodes exist (they are conditioned averages, as the case 0/0 is not defined).

**w_s_io.hist and w_k_io.hist files**

This files contain the (one realization, not averaged over the ensemble) graph-average existing occupation number as function of the product of strenghts or degrees in the format:
```
x^{out}_i*x^{in}_j \bar{t_ij} \std{\bar{t_ij}}
```
Where x stands for degree and strengths respectively and the averages $\bar{x}$ are taken over a single graph realization (and the standard deviations also).

**_w.hist files**
This files contain the distribution of existing occupation numbers in the format:
```
bin_id bin_min bin_max Bin_count Bin_std CCDF
```
For the case of a single realization, the histogram is not normalized, while for the ensemble average "ens_" file the histogram is normalized.

**Example**:

For the out-strength distribution, for instance, 3 files are produced for 10 reps and *avs=500*:

*N5000avs499.88180_w.hist*			One realization of the network P(w)
*N5000avs499.88180_ens_r1_w.hist*	Ensemble average distribution <P(w)>
*N5000avs499.88180_ens_r1node.list*	Ensemble averages for different node features
*N5000avs499.88180node.list*		Graph values for single realization over node features
*N5000avs499.88180_w_s_io.hist*		Existing Occupation number average as function of s_in s_out over single realization
*N5000avs499.88180_w_k_io.hist*		Existing Occupation number average as function of k_in k_out over single realization

The undirected version of the algorithm produces equivalent files for the *in* and *out* values.


## Acknowledgements

We would like to thank Pol Colomer and Sergio Oller for their very useful comments and suggestions.



## License

Copyright 2013 Oleguer Sagarra and Francesc Font-Clos.
All rights reserved. 
Code under License GPLv3.

## Roadmap

The software is for the moment ready and works smoothely under apropiate conditions.
It could be extended to include all the cases described in ref. 1, but this is work in progress,
If you wish to contribute please contact osagarra@ub.edu.



## Known issues 

Don't be too brave simulating high average occupation numbers if N_nodes is high, you'll probably reach the limit of your computer (specially with the Canonical Ensemble). 
The algorithm has been tested up to 5e5 nodes but not further.

## Additional notes

Notes:

Check the DOCS/requirements to see the dependencies

This soft has been tested on Linux Ubuntu 12.10 and MacOS X 6.x

Recommended compiler for mac is Clang while icc for Linux.


