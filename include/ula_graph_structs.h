#ifndef ULA_GRAPH_STRUCTS_H
#define ULA_GRAPH_STRUCTS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

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

#endif
