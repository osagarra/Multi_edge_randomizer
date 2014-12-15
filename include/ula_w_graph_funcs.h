#ifndef ULA_W_GRAPH_FUNCS_H
#define ULA_W_GRAPH_FUNCS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_fit.h>
#include "ula_gen_funcs.h"
#include "ula_stat_funcs.h"
#include "ula_graph_structs.h"


/******* Allocation ********/
w_graph* w_graph_alloc(int N_nodes);
void w_graph_free(w_graph* node, int N_nodes);

/******* Add edge ********/
void w_graph_add_multi_link(w_graph * node, int N_nodes, int origin, int dest, int weight);
void w_graph_add_multi_link_undirected(w_graph * node, int N_nodes, int origin, int dest, int weight);

/******* Printing ********/
void w_graph_print_adj_list(w_graph* node, int N_nodes, char* output);

/******* S & k's ********/
int ** w_graph_compute_s(w_graph* node, int N_nodes);
int ** w_graph_compute_k(w_graph* node, int N_nodes);
double ** w_graph_compute_k_analitic(w_graph* node, int N_nodes, int self_opt);
double ** w_graph_compute_k_analitic_from_s_directed(int** s, int N_nodes, int self_opt);
double ** w_graph_compute_k_analitic_from_s_undirected(int* s, int N_nodes, int self_opt);

/******* Snn & knn's ********/
double ** w_graph_compute_s_nn(w_graph* node, int N_nodes, int weight, int opt_dir);
/******* Y2 ********/
double ** w_graph_compute_Y2(w_graph * node, int N_nodes, int opt_dir);
//double ** w_graph_compute_xy(w_graph * node, int N_nodes);

/******* Clsutering ********/
double ** w_graph_compute_clust(w_graph * node, int N_nodes);
/******* w's ********/
int * w_graph_compute_w(w_graph* node, int N_nodes, int* aux, int zeros);
double** w_graph_compute_p_w_analitic_from_s_undirected(int maxt, double binn, int* s, int N_nodes, int self_opt, int* len);
double** w_graph_compute_p_w_analitic_from_s_directed(int maxt, double binn, int** s, int N_nodes, int self_opt, int* len);
double * w_graph_compute_w_ss(w_graph* node, int N_nodes, int weight);
/******* Distances *****/
//double ** w_graph_compute_dists(w_graph* node, int N_nodes);
/*******  Entropies *****/
double w_graph_entropy(w_graph* node, int N_nodes);
void w_graph_print_entropy(double* seq,  int len,char* output);
double w_graph_loglikelyhood(w_graph* node,int N_nodes, double** wij);
/******* All_stats *****/
void w_graph_node_stats_list(w_graph* node, int N_nodes, int run, double av_k, int opt_dir, int opt_clust, int self_opt);
int w_graph_total_weight( w_graph* node, int N_nodes);
int w_graph_total_edges( w_graph* node, int N_nodes);
void w_graph_all_stats(w_graph* node, int N_nodes, int run, double bin_exp,double av_k, int opt_dir, int self_opt, int w_anal);

/******* Ensemble stats *****/
void w_graph_node_stats_ensemble(w_graph* node, int N_nodes, double** container, double** container2, int** node_nonzero, double* T_container, int opt_dir, int opt_clust );
void w_graph_node_stats_ensemble_print(int reps, int N_nodes, double* Tcont, double** cont, double ** cont2, int** node_nonzero, double av_k, double bin_exp, int len_acc, int opt_dir);

gsl_histogram ** w_graph_all_stats_ensemble_allocate(int dir, int s_min, int s_max, int k_min, int k_max, int w_max);
void w_graph_all_stats_ensemble_update(gsl_histogram** acc, w_graph* node, int N_nodes, int dir);
void w_graph_all_stats_ensemble_print(gsl_histogram** acc, int len, int reps, int N_nodes, double av_k, int opt_dir);
/****************************************************************************
 * Graph algebra *
 ****************************************************************************/
void w_graph_sum_graphs(w_graph* node1, w_graph* node2, int N_nodes);
void w_graph_normalize(w_graph* node, int reps, int N_nodes);

#endif