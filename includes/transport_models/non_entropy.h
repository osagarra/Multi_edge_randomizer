#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>

#include "../ula_funcs.h"
#include "../graph_funcs.h"
#include "../w_graph_funcs.h"
#include "../../main_simus.h"

double compute_gamma_frenchies(double surface);
double compute_p_french(int *s_in, double **d, double gamma, int dest, int origin, int* ind_in , int * ind_out, int len, int self_opt, int verbose);
void check_norm_french(int *s_in, double **d, double gamma, int origin, int* ind_in, int *ind_out, int len, int self_opt, int verbose);
w_graph* w_graph_seq_gravity_directed_graph(int N_nodes,int **s, double **d, int Reps, double gamma, gsl_rng* randgsl, int Max_fails, int self_opt, int verbose);
/***********/
int destination_rad_model(int** s, double** d, int origin, int N_nodes, int T, gsl_rng * randgsl);
void radiation_model(int **w_real,int N_nodes,int **s, double **d, int Reps,int seed, char *output_name, char *output_indices);
void radiation_model_multinomial(int **w_real,int N_nodes,int **s, double **d, int Reps,int seed, char *output_name, char *output_indices);
int compute_sij_rad(int  **s, double **dist, int origin, int dest, int N_nodes);
/***********/
void bose_model_dist(int **w_real,int N_nodes,int **s, double **x, double **d, double gamma, int Reps,int seed, char *output_name, char *output_indices);
void bose_model(int **w_real,int N_nodes,int **s, double **x, int Reps,int seed, char *output_name, char *output_indices);
/***********/
void wilson_model_distance(int **w_real,int N_nodes,int **s, double **x, double **d, double gamma, int Reps,int seed, char *output_name, char *output_indices);
void wilson_model(int **w_real,int N_nodes,int **s, double **x,int Reps,int seed, char *output_name, char *output_indices);
