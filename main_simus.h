/********************************************************************************
 * 
 *                         Frenchies all
 * 
 * ******************************************************************************/

#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cblas.h>
#include <unistd.h>


typedef struct w_graph{
	int idnum;
	int kin;
	int kout;
	int sin;
	int sout;
	double x;
	double y;
	double loc_x;
	double loc_y;
	int *out;
	int *in;
	int mem_in;
	int mem_out;
	int* w_out;
	int* w_in;
}w_graph;




#include "./includes/ula_funcs.h"
//#include "./includes/solver.h"
#include "./includes/stat_funcs.h"
#include "./includes/graph_funcs.h"
#include "./includes/w_graph_funcs.h"
//#include "./includes/transport_models/transport_models.h"
#include "./includes/transport_models/null_models.h"




#endif
