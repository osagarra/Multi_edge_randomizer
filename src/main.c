/********************************************************************************
 * 
 *			GENERATE Maximum entropy ensembles for networks with given <s> = T/(N-1)
 *
 *	This program generates networks (and statistics) with pre-defined strength distribution using either a
 *  macro-canonical approach (poisson and multinomial) or a canonical approach or Micro-canonical approach.
 *
 * 	See README for details and the DOCS folder
 *
 * 
 * Author: Oleguer Sagarra and Francesc Font Clos, 2013.
 * 
 * 
 * Check below for information on the program or type ./MultiEdgeGen without arguments for help
 *
 *	Output:
 *		
 *		*.tr: A file stored in list format of flows between given nodes
 *			node i node j t_ij
 *		*.hist: Diverse statistics on the network (indicative names on the files)
 *		*.list: Node attributes average with average deviation (see files).
 * 
 * ******************************************************************************/




#include "main.h"


int main(int argc, char *argv[]){
	printf(\
	"################################\n"
	"########## Maximum entropy multi-edge networks ###########\n"
	"####################################\n");

/***********************************************************************
 default params
************************************************************************/
	//int  	N_nodes				= atoi(argv[1]);
	//int  	Reps				= atoi(argv[2]);
	int  	seed      			= 1;
	int opt_dir=-1;
	int N_nodes=-1;
	//int opt_dir					= atoi(argv[3]);
	int meth					=2;
	int print_tr				=0;
	double bin_exp				=-1;
	int reps					=100;
	int verbose					=0;
	int opt_clust				=0;
	int self_opt				=1;
	int w_anal                  =0;
	int header 					=1;
	char* file_s;
	int ch;
/***********************************************************************
 read args
************************************************************************/
	        while ((ch = getopt(argc, argv, "N:s:d:e:p:f:x:r:v:c:l:w:h:")) != -1) {
	             switch (ch) {
	             case 's': /* seed */
	                     seed=atoi(optarg);
	                     break;
			     case 'N': /* N_nodes */
			             N_nodes=atoi(optarg);
			             break;
	             case 'd': /* dir_opt */
	                     opt_dir=atoi(optarg);
	                     break;
	             case 'e': /* ensemble opt */
	                     meth=atoi(optarg);
	                     break;
	             case 'p': /* print opt */
	                     print_tr=atoi(optarg);
	                     break;
	             case 'f': /* file strength */
	                     file_s=optarg;
	                     break;
	             case 'x': /* exp */
	                     bin_exp=atoi(optarg);
	                     break;
			     case 'r': /* reps */
	                     reps=atoi(optarg);
			             break;
			     case 'v': /* verbose */
			             verbose=atoi(optarg);
			             break;
			     case 'c': /* clust */
			             opt_clust=atoi(optarg);
			             break;
			     case 'l': /* self loop */
			             self_opt=atoi(optarg);
			             break;
				 case 'w': /* w_analitical */
					 	 w_anal = atoi(optarg);
					     break;
				 case 'h': /* header */
				 		 header = atoi(optarg);
						 break;
	             default:
				 {
	                     printf("Unknown flag %c\n", ch);
	                     exit(EXIT_FAILURE);
	             }
			 	}
	}
	if((N_nodes<0)||(opt_dir<0))
	{
 		fprintf(stderr,	"\nCorrect usage is: ./simus -args \n\nWhere:\n\n"
 				" *  Compulsory items:\n"
 				" *        -N N_nodes. Number of nodes (int)\n"
 				" *        -d dir_opt. Undirected (0) or Directed (1)\n"
 			    " *        -f file_s Path to file with strength sequence in form on each line: node_num(int) s_out(int) s_in(int) in the directed case, node_num (int) s(int) otherwise\n"
 				" *  Optional items: \n"				
 				" *        -s seed.initial seed for random generator (int) (default=1)\n"
 				" *        -e ensemble_opt. Method: 0 (canonical, multinomial), 1 (grand-canonical, poisson + multinomial), 2 (grand-canonical, poisson indep.), 3(micro-canonical)\n"
 				" *        (Default=2) \n"
 				" *        -p print_opt. Print adj list? 0 (no), 1 (yes) (Default=0)\n"
 				" *        -x Exponent for log-binning (-1 for no log binning) (Default=-1)\n"
 				" *        -r Number of reps for averaging (int) (Default=100)\n"
 				" *        -v Verbose (1 for on, 0 for off) (Default 0)\n"
 				" *        -c Clustering option (1 for yes) (warning: Depending on av_s makes simulations orders of magnitude slower) (Default=0)\n"
 				" *        -l Self-loop option (>0 for accepting them) (Default =1) \n"
                " *        -w Compute analytic distribution of weights? (>0 for yes, takes some time) (Default=0)\n\n"
				" *        -h Number of header lines on file_s (default=1)\n\n"
				"Please, read the DOCS/README file for more info!\n");
 		return 0;
 	}

	/****** Check all in params are good ******/
	if(bin_exp<=1) bin_exp=1.05;
    if((opt_dir!=1)&&(opt_dir!=0))
    {
		printf("Select directed or undirected!, aborting...\n");
		abort();
	}
	if((meth>3) || (meth<0))
	{
		printf("Select apropiate method (canonical, grandcanonical...), aborting ....\n");
		abort();
	}
/*
	if((opt_dir==1) && (meth==2))
	{
		printf("Sorry, not implemented the grand canonical for poisson independent directed case, aborting...\n");
		abort();
	} 
*/
	if((print_tr != 0) && (print_tr != 1))
	{
		printf("Select if you want to print the adj list, aborting ....\n");
		abort();
	}
	if(opt_clust==1)
	{
		if(opt_dir==1)
		{
			printf("Ignoring clustering option, only defined for undirected networks ....\n");
			opt_clust=0;
		}else{
			printf("CLustering option selected! This may cause low performance for high <s>! ....\n");
		}
	}
  /*** Set rand generator (uses GLS THAU) ***/ 
	if(seed==1)
	{
		seed = seed + time(NULL); // Change for trully random generation
		if(verbose>1) printf("Default seed given, using processor time...\n");
	}
	gsl_rng * randgsl = gsl_rng_alloc(gsl_rng_taus);	/// we initialize the random generator of the gsl
	gsl_rng_set(randgsl,seed);

 /*** Print some info ***/

	printf("SEED=%i\n",seed);
	printf("Average over %d reps\n", reps);
	//printf("Iterations=%i\n",Reps);


	int T;
	//long double Tfake;
	double ** x2;
	double * x;
	int** xx2;
	int* xx;
	double av_k;
 /***********************************************************************
 	Allocating memory + reading distro
 ************************************************************************/ 	
	if(opt_dir==1)
	{
		xx2 = read_node_list_int(file_s, N_nodes,header); // strenght sequence (ints)
		T=sum_vec_int(xx2[0],N_nodes); // T is \sum_i s_i (for all cases)
		if(meth!=3)
		{
			x2 = cast_mat_double(2,N_nodes); // doubles
			x2[0] = vec_int_to_double(xx2[0],N_nodes);
			x2[1] = vec_int_to_double(xx2[1],N_nodes);
		}
	}else{
		xx = read_node_list_int_undir(file_s, N_nodes,header); // strenght sequence (ints)
		T=sum_vec_int(xx,N_nodes); // T is \sum_i s_i (for all cases)
		if(T%2!=0)
		{
			/*T++;
			if((seed>=0) && (seed<N_nodes))
			{
				xx[seed]++;
			}else{
				xx[0]++;
			}
            
			if(verbose>0)
			{
				//printf("Warning: T is not even, added an additional stub at seed node!\n");
                printf("Warning: T is not even, aborting\n");
                abort();
			}
			*/
		}
		if(meth!=3)
		{
			x = vec_int_to_double(xx, N_nodes); // convert to doubles
		}
	}
	av_k=(double)T/(double)N_nodes;
	//Tfake = (long int)T;
	if(verbose==1)
	{
		if(opt_dir==1)
		{
			printf("- T to be sorted: %d -> <s>=%f || s_max/T : (out) =%f (in) %f\n",T,(double)T/(double)N_nodes, (double)max_value_int(xx2[0],N_nodes)/(double)T, (double)max_value_int(xx2[1],N_nodes)/(double)T );
		}else{
			printf("- T to be sorted: %d -> <s>=%f || s_max/2T = %f\n",T,(double)T/(double)N_nodes, (double)max_value_int(xx,N_nodes)/(double)(T) );
		}
		if(self_opt>0)
		{
			printf("Warning! Self-loops accepted! \n");
			fflush(stdout);
		}
	}
	double * ps;
/***********************************************************************
	Preparing ensemble reps!
 ************************************************************************/ 	
	int r, len_acc_nodes;
	W_GRAPH* WG;
	double ** node_cont; //node info container for averaging
	double ** node_cont2; // node info std for averaging
	int ** node_nonzero; // node count non-zero strength/degree
	double * Tcont; // Total T container (av and std)
	int s_max, w_max;
	//int s_min;
	if(opt_dir==1)
	{
		len_acc_nodes=14; // in-out variables
		s_max=(int)max_value_int(xx2[0],N_nodes)*5;
		//s_min=min_value_int(xx2[1],N_nodes)/3.;
	}else{
		if(opt_clust==1)
		{
			len_acc_nodes=9;
		}else{
			len_acc_nodes=7;
		}
		s_max=(int)max_value_int(xx,N_nodes)*5;
		//s_min=(double)min_value_int(xx,N_nodes)/3.;
	}
	if((int)s_max<1)
	{
		s_max=10;
	}
	w_max=(int)(((double)s_max)*((double)s_max/(double)T));
	if(w_max<1)
	{
		w_max=2;
	}
	node_cont=cast_mat_double(N_nodes,len_acc_nodes);
	node_cont2=cast_mat_double(N_nodes,len_acc_nodes);
	node_nonzero=cast_mat_int(N_nodes,2); // to count the number of times a node is non-connected (in-out)
	Tcont=cast_vec_double(2); // mean and std deviation
	gsl_histogram ** acc_ensemble= w_graph_all_stats_ensemble_allocate(opt_dir, 0, 0, 0, 0, w_max);
	

	/***********************************************************************
		Fill analytic calculation for degree (and optionally clustering)
	 ************************************************************************/ 	

	int i;
	if(opt_dir==1)
	{
		double** k=w_graph_compute_k_analitic_from_s_directed(xx2,N_nodes, self_opt);
		for(i=0;i<N_nodes;i++)
		{
			node_cont[i][1]=k[0][i]*reps;
			node_cont2[i][1]=(k[2][i]+k[0][i]*k[0][i])*reps; // directly sigma!
			node_cont[i][7]=k[1][i]*reps;
			node_cont2[i][7]=(k[3][i]+k[1][i]*k[1][i])*reps; // directly sigma!
		}
		free_mat_double(k,4);
	}else{
		double ** k=w_graph_compute_k_analitic_from_s_undirected(xx,N_nodes, self_opt);
		for(i=0;i<N_nodes;i++)
		{
			node_cont[i][1]=k[0][i]*reps;
			node_cont2[i][1]=(k[1][i]+k[0][i]*k[0][i])*reps; // snd moment
		}
		free_mat_double(k,2);
	}

	
	/***********************************************************************
		Preparing the distribution of occupation numbers
	 ************************************************************************/ 	
		/** Probabilities choice (could be extended to other types) **/
	if (meth<3)
	{
		if(opt_dir==1)
		{
			if (meth==2) 
			{
				free_mat_int(xx2,2); // saves memory
			}else{
				ps = prob_mult_s_dir(x2, N_nodes, self_opt);
			}
		}else{
			if (meth==2)
			{
				free(xx); // saves memory
			}else{
				ps = prob_mult_s_undir(x, N_nodes, self_opt);
			}
		}
	}
/***********************************************************************
	Start of ensemble reps!
 ************************************************************************/

	for(r=0;r<reps;r++)
	{
		if(meth==2)
		{
			if(verbose==1) printf("============### Poisson model ####===========\n"); fflush(stdout);
			if(opt_dir==1)
			{
				//printf("-- Sorry, not implemented poisson for grand-canonical directed case, aborting --\n");
				//abort();
				WG = uncorrelated_poisson_directed_graph2(x2, N_nodes, randgsl, verbose, self_opt);
				// not implemented!
			}else{
				WG = uncorrelated_poisson_undirected_graph2(x, N_nodes, randgsl, verbose, self_opt); // more stable with big T,N

			}
		}else if(meth==3){
			if(verbose==1) printf("============### Rewiring model ####===========\n"); fflush(stdout);
			if(opt_dir==1)
			{
				WG = uncorrelated_computational_directed_graph(xx2,N_nodes, randgsl,verbose, self_opt);
			}else{
				WG = uncorrelated_computational_undirected_graph(xx,N_nodes, randgsl, 100000, verbose, self_opt);
			}
		}else{
			if(meth==0)
			{
				if(verbose==1)printf("============### Multinomial model ####===========\n"); fflush(stdout);
				if(opt_dir==1)
				{
					WG = uncorrelated_multinomial_directed_graph(ps, x2, N_nodes, T, randgsl,verbose, self_opt);
				}else{
					WG = uncorrelated_multinomial_undirected_graph(ps, x, N_nodes, T, randgsl,verbose, self_opt);
				}
			}else{
				if(verbose==1)printf("============### Poisson Multinomial model ####===========\n"); fflush(stdout);
				if(opt_dir==1)
				{
					WG = uncorrelated_poisson_multinomial_directed_graph(ps, x2, N_nodes, T, randgsl, verbose, self_opt);
				}else{
					WG = uncorrelated_poisson_multinomial_undirected_graph(ps, x, N_nodes, T, randgsl, verbose, self_opt);
				}
			}
		}
		if(r==0) // if first rep, store all stats
		{
			w_graph_all_stats(WG, N_nodes, 0,  bin_exp, av_k, opt_dir,self_opt, w_anal);
			w_graph_node_stats_list(WG,N_nodes,0, av_k, opt_dir, opt_clust, self_opt);
            char cadena[100];
            if(print_tr==1)
            {
                if(verbose>0)
                {
                    printf("Printing adj matrix\n");
                    fflush(stdout);
                }
                sprintf(cadena,"N%d_cust.tr",N_nodes);
                w_graph_print_adj_list(WG, N_nodes, cadena);
            }
		}
		w_graph_node_stats_ensemble(WG,N_nodes,node_cont,node_cont2,node_nonzero, Tcont,opt_dir,opt_clust);
		w_graph_all_stats_ensemble_update(acc_ensemble,WG, N_nodes, opt_dir);
		w_graph_free_nodes(WG, N_nodes);
		free(WG);
	}
	if(meth<2) free(ps);
/***********************************************************************
	 Print output
************************************************************************/ 	
	int len=2;
	/*
	if(opt_dir==1)
	{
		len=10;
	}else{
		len=6;
	}
	*/
	w_graph_all_stats_ensemble_print(acc_ensemble, len, reps, N_nodes,av_k, opt_dir);
	w_graph_node_stats_ensemble_print(reps, N_nodes, Tcont, node_cont, node_cont2, node_nonzero, av_k, bin_exp,len_acc_nodes, opt_dir);
	gsl_rng_free (randgsl);
	return 0;
}
