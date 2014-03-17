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

[2]

This software is distributed as it is and may be used for scientific uses only,
please report bugs to osagarra@ub.edu.

If you do use our software, please cite us!


############# Install ################

- Just unpack the .tar package and run:
	./compile_simus.sh for compilation (./compile_simus_mac.sh for mac users with Clang compilers).

Notes:

- Check the DOCS/requirements to see the dependencies

- This soft has been tested on Linux Ubuntu 12.10 and MacOS X 6.x

- Recommended compiler for mac is Clang while icc for Linux.

########### Usage ###########

1. To provide yourself the strength sequence and randomize the network using a maximum entropy approach, type:

./simus_cust 
*	Arguments:
 *		1. Number of nodes (int)
 *		2.initial seed for random generator (int)
 * 		3. Undirected (0) or Directed (1)
 * 		4. Method: 0 (canonical, multinomial), 1 (grand-canonical, poisson + multinomial), 2 (grand-canonical, poisson indep.), 3(micro-canonical)
 * 		5. Print adj list? 0 (no), 1 (yes)
 * 		6. file_s --> Path to file with strength sequence in form on each line: "node_num(int) s_out(int) s_in(int)\n" in the directed case, "node_num (int) s(int)\n" otherwise
 *		7. Exponent for log-binning (-1 for no log binning)
 *		8. Number of reps for averaging (int)
 *		9. Verbose (1 for on, 0 for off)
 * 		10. Clustering option (warning: Depending on av_s makes simulations orders of magnitude slower)
*		11. Self-loop option (1 for accepting self-loops, 0 for non accepting)


For both cases:

The generated outputs are files with extension "*.hist" if they correspond to network statistics (log binned), "*.tr" for the network in weighted adjacency format and "*.list" for the node list with different features.

For the out-strength distribution, for instance, 3 files are produced for exp=2.5, 10 reps and avs=500.:

N5000avs499.88180expo-1.00_sOUT.hist # One realization of the networks
N5000avs499.88180expo-1.00_ens_r1_sOUT.hist # Ensemble average distribution <P(s)>
N5000avs499.88180expo-1.00_ens_r1node.list # Ensemble averages for different node features

The undirected version of the algorithm produces equivalent files for the "in" and "out" *.hist

The file *ens_rXnode_list.list contains ensemble averages for different node features, together with their standard deviations.
# Node_num  k sigma_k   k_analitic sigma_k_analitic  s sigma_s   Y2 sigma_Y2   k_nn sigma_k_nn   k^w_nn sigma_k^w_nn  s^w_nn sigma_s^w_nn (optionally  c sigma_c  c^w sigma_c_w) #

Please note that the averages for Y2, Clustering and Average neighbor properties are performed only over realizations where the nodes exist!


######## Roadmap ##########

The software is for the moment ready and works smoothely under apropiate conditions.

It could be extended to include all the cases described in ref. 1, but this is work in progress,

If you wish to contribute please contact osagarra@ub.edu.



####### Known issues #######

Don't be too brave simulating high average occupation numbers if N_nodes is high, you'll probably reach the limit of your computer. The algorithm has been tested up to 5e5 nodes and av_w = 0.1 but not further.
