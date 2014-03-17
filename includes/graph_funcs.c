/************************************************************
 *
 *                    Graph Library
 *
 *		Functions useful when dealing with spatial directed weighted networks
 *
 *
 *************************************************************/



#include "graph_funcs.h"

/****************************************************************************
 *  Reading functions *
 ****************************************************************************/

int** read_edge_list(char *input_name, int num_nodes){
	printf("reading edge list file and converting to Weighted adjacency matrix...\n");
	FILE* input=open_file("r",input_name);
	int** d=cast_mat_int(num_nodes,num_nodes);
	int i,j,dij;
	int n=0;
	while (!feof(input))
	{		///we start reading	
		if(fscanf(input, "%d %d %d\n", &i, &j, &dij)!=3)
		{
			printf("error al llegir\n");
			abort();
		}
  		d[i][j]=dij;
  		n+=dij;
  		//printf("%i",n);fflush(stdout);
  	}
	printf("Total num of trips %i\n",n);
	fclose(input);
	return d;
}

double*** read_edge_list_double(char *input_name, int num_nodes){ // check
	printf("reading average edge list file and converting to Weighted adjacency matrix...\n");
	FILE* input=open_file("r",input_name);
	double** d=cast_mat_double(num_nodes,num_nodes);
	double** d2=cast_mat_double(num_nodes,num_nodes);
	int i,j;
	double dij,dij2;
	double n=0;
	while (!feof(input))
	{		///we start reading	
		if(fscanf(input, "%d %d %lf %lf\n", &i, &j, &dij, &dij2)!=4)
		{
			printf("error al llegir\n");
			abort();
		}
  		d[i][j]=dij;
  		d2[i][j]=dij2;
  		n+=dij;
  		//printf("%i",n);fflush(stdout);
  	}
	printf("Total average num of trips %lf\n",n);
	fclose(input);
	double *** dtot=(double***)malloc(sizeof(double**)*2);
	dtot[0]=d;
	dtot[1]=d2;
	return dtot;
}

/************************************/

int** read_node_list_int(char *input_name,int num_nodes){
	printf("reading node attribute file...\n");
	FILE* input=open_file("r",input_name);
	int** s=(int**)malloc(sizeof(int*)*2);
	s[0]=(int*)malloc(sizeof(int)*num_nodes);
	s[1]=(int*)malloc(sizeof(int)*num_nodes);
	int n,out,in;
	int m=0;
	int att_in,att_out;
	int flag=0;
	att_in=att_out=0;
	while (!feof(input)){		///we start reading	
		if(fscanf(input, "%d %d %d\n", &n, &out, &in)!=3)abort();
		//assert(fscanf(input, "%d %d %d\n", &n, &out, &in)==3);
  		if(m!=n)flag=1;
  		s[0][m]=out;
  		s[1][m]=in;
  		m++;
  		att_in+=in;
  		att_out+=out;
 		}
	if(num_nodes!=m){
		printf("Check your node list, its different (size %d) than the declared number of nodes (%d)!\n",m,num_nodes);
		abort();
	}else if(flag==1){printf("Node numbers are not congruent with given list!");
	}
	printf("Total number of nodes %i \t Total in_att %i \t Total out_att %i Difference %i\n",m,att_in,att_out,abs(att_in-att_out));
	fclose(input);
	return s;
}

int* read_node_list_int_undir(char *input_name,int num_nodes){
	printf("reading node attribute file...\n");
	FILE* input=open_file("r",input_name);
	int* s=(int*)malloc(num_nodes*sizeof(int));
	int n,out;
	int m=0;
	int att_out;
	int flag=0;
	att_out=0;
	while (!feof(input)){		///we start reading	
		if(fscanf(input, "%d %d\n", &n, &out)!=2)abort();
		//assert(fscanf(input, "%d %d %d\n", &n, &out, &in)==3);
  		if(m!=n)flag=1;
  		s[m]=out;
  		m++;
  		att_out+=out;
 		}
	if(num_nodes!=m){
		printf("Check your node list, its different (size %d) than the declared number of nodes (%d)!\n",m,num_nodes);
		abort();
	}else if(flag==1){printf("Node numbers are not congruent with given list!");
	}
	printf("Total number of nodes %i \t Total att %i \n",m,att_out);
	fclose(input);
	return s;
}



/************************************/

double** read_node_list_double(char *input_name,int num_nodes){
	printf("reading node attribute file...\n");
	FILE* input=open_file("r",input_name);
	double ** s=(double**)malloc(sizeof(double*)*2);
	s[0]=(double*)malloc(sizeof(double)*num_nodes);
	s[1]=(double*)malloc(sizeof(double)*num_nodes);
	int n;
	int m=0;
	double att_in,att_out,in,out;
	int flag=0;
	att_in=att_out=0;
	while (!feof(input)){		///we start reading	
		if(fscanf(input, "%d %lf %lf\n", &n, &out, &in)!=3)abort();
		//assert(fscanf(input, "%d %lf %lf\n", &n, &out, &in)==3);
  		if(m!=n)flag=1;
  		s[0][m]=out;
  		s[1][m]=in;
  		m++;
  		att_in+=in;
  		att_out+=out;
 		}
	if(num_nodes!=m){
		printf("Check your node list, its different (size %d) than the declared number of nodes (%d)!\n",m,num_nodes);
		abort();
	}else if(flag==1){printf("Node numbers are not congruent with given list!");
	}
	printf("Total number of nodes %i \t Total in_att %lf \t Total out_att %lf Difference %lf\n",m,att_in,att_out,fabs(att_in-att_out));
	fclose(input);
	return s;
}




