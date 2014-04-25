############################################################
##########   Maximum entropy multi-edge network   ##########
##########              generator                 ##########
##########                                        ##########
############################################################

Authors: Oleguer Sagarra and Francesc Font Clos
Year: 2013
Details on the theory (please cite us if you use the software):

[1] Statistical mechanics of multiedge networks, Sagarra O., Pérez-Vicente C. and Díaz-Guilera, A.  Phys. Rev. E 88, 062806 (2013)

    http://pre.aps.org/abstract/PRE/v88/i6/e062806

[2] The effect of fixing node strengths on multi-edge networks, Sagarra O., Font-Clos F., Pérez-Vicente C. and Díaz-Guilera, A. Arxiv pre-print.

	http://arxiv.org/abs/1404.3697

This software is distributed as it is and may be used for scientific uses only,
please report bugs to osagarra@ub.edu.

If you do use our software, please cite us!

Additional details on the implementation can be found in the DOCS/details.pdf

############# Install ################

- Just unpack the .tar package and run:
	./compile_simus.sh for compilation (./compile_simus_mac.sh for mac users with Clang compilers).



Notes:

- Check the DOCS/requirements to see the dependencies

- This soft has been tested on Linux Ubuntu 12.10 and MacOS X 6.x

- Recommended compiler for mac is Clang while icc for Linux.

########### Usage ###########

1. To provide yourself the strength sequence and randomize the network using a maximum entropy approach, type:

./simus_cust help to list all the available options

There are options for accept/not self-loops, choose the ensemble to use (canonical, grand canonical, micro-canonical) and directed/undirected mode. Other options are also available.


########### Inputs ##########

The input needed is a strength sequence. Use the flag -h to specify the number of lines to skip while reading the file (header files).

For the undirected case:
	The format needs to be "%d %d\n" where the first index is the node number, the second element being the strength
For the directed case:
	The format needs to be "%d %d %d\n" where the first index is the node number, the second element being the strength_out and the last being the strength_in.

############ Outputs #########

For both directed and undirected cases:

The generated outputs are files with extension "*.hist" if they correspond to network statistics (log binned), "*.tr" for the network in weighted adjacency format and "*.list" for the node list with different features.


The files with the suffix "ens" reffer to averages over ensembles. The ones without the "ens" suffix refer to a single realization (and hence the .list files do not contain averages nor std, but only raw values for each node).


- .list files:
   The file *ens_rXnode_list.list contains ensemble averages for different node features, together with their standard deviations.
   For the undirected case: (Note that the clustering only appears if -c flag is set)	
# Node_num  k k_std k_analitic k_analitic_std  s s_std   Y2 Y2_std   k_nn k_nn_std   k^w_nn k^w_nn_std  s^w_nn s^w_nn_std  (optionally  clust c_std  c^w c^w_std)

Explicit formulas:
	k_i = \sum \Theta(t_ij), k^{anal}_i = \sum (1-e^{-s_is_j/T}), s_i = \sum t_{ij}, Y2 = \sum t_ij^2 / s_i
	k_nn_i = \sum \Theta(t_{ij}) k_j / k_i, k^w_nn_i = \sum t_{ij} k_j /s_i, s^w_nn = \sum t_{ij} s_j / s_i
	c = \sum_{jk} \Theta(t_ij)\Thetat(t_jk)\Theta(t_ki) / k_i(k_i-1)
	c^w = \sum_{jk} (t_ij + t_ik) \Theta(t_ij) \Thetat(t_jk)\Theta(t_ki) / [2s_i(k_i-1)]

For the directed case: (In this case the clustering is not defined)

# Node_num  k k_std k_analitic k_analitic_std  s s_std   Y2 Y2_std   k_nn k_nn_std   k^w_nn k^w_nn_std  s^w_nn s^w_nn_std 

First for the out case, then for the in case. Explicit formulas:
Out:
	k_i^{out} = \sum_j \Theta(t_ij), k^{anal}_i^{out} = \sum_j (1-e^{-s_i^out}s_j^{in}/T}), s_i^{out} = \sum_j t_{ij}, 
	Y2^{out} = \sum_j t_ij^2 / s_i^{out}
	k_nn_i^{out} = \sum_j \Theta(t_{ij}) k_j^{in} / k_i^{out}, k^w_nn_i = \sum_j t_{ij} k_j^{in} /s_i^{out}, 
	s^w_nn = \sum_j t_{ij} s_j^{in} / s_i^{out}
In:
	k_j^{in} = \sum_j \Theta(t_ij), k^{anal}_j = \sum_j (1-e^{-s_i^out s_j^in/T}), s_j^in = \sum_i t_{ij}, 
	Y2^{in} = \sum_i t_ij^2 / s_j^{in}
	k_nn^{in} = \sum_i \Theta(t_{ij}) k_i^{out} / k_j^{in}, k^w_nn = \sum_i t_{ij} k_i^{out} /s_j^{in}, 
	s^w_nn^{in} = \sum_i t_{ij} s_i^{out} / s_j^{in}
	

Please note that the averages for Y2, Clustering and Average neighbor properties are performed only over realizations where the nodes exist!

- w_s_io.hist and w_k_io.hist files:

This files contain the (one realization) average existing occupation number as function of the product of strenghts or degrees in the format:
# x^{out}_i*x^{in}_j \bar{t_ij} \std{\bar{t_ij}}
Where x stands for degree and strengths respectively and the averages $\bar{x}$ are taken over a single graph realization (and the standard deviations also).

- Files _w.hist files:

This files contain the distribution of existing occupation numbers in the format:
# bin_id bin_min bin_max Bin_count Bin_std CCDF

For the case of a single realization, the histogram is not normalized, while for the ensemble average "ens_" file the histogram is normalized.


--- For example: ---

For the out-strength distribution, for instance, 3 files are produced for 10 reps and avs=500.:

N5000avs499.88180_w.hist # One realization of the network P(w)
N5000avs499.88180_ens_r1_w.hist # Ensemble average distribution <P(w)>
N5000avs499.88180_ens_r1node.list # Ensemble averages for different node features
N5000avs499.88180node.list # Graph values for single realization over node features
N5000avs499.88180_w_s_io.hist # Existing Occupation number average as function of s_in s_out over single realization
N5000avs499.88180_w_k_io.hist # Existing Occupation number average as function of k_in k_out over single realization


The undirected version of the algorithm produces equivalent files for the "in" and "out" *.hist


######## Roadmap ##########

The software is for the moment ready and works smoothely under apropiate conditions.

It could be extended to include all the cases described in ref. 1, but this is work in progress,

If you wish to contribute please contact osagarra@ub.edu.



####### Known issues #######

Don't be too brave simulating high average occupation numbers if N_nodes is high, you'll probably reach the limit of your computer. The algorithm has been tested up to 5e5 nodes and av_w = 0.1 but not further.
