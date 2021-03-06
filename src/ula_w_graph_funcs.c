/************************************************************
 *
 *                    w_Graph Library
 *
 *		Functions useful when dealing with directed weighted networks in edge list format
 *
 *
 *************************************************************/



#include "ula_w_graph_funcs.h"

/****************************************************************************
 * Allocation *
 ****************************************************************************/
W_GRAPH* w_graph_alloc(int N_nodes){
    int i;
    W_GRAPH* WG = (W_GRAPH*)malloc(sizeof(W_GRAPH));
    WG->node = (NODE*)malloc(sizeof(NODE)*N_nodes);
    WG->N_nodes = N_nodes;
    for(i=0;i<N_nodes;i++)
    {
        WG->node[i].idnum=i;
        //node[i].x=x[0][i];
        //node[i].y=x[1][i];
        WG->node[i].kin=0;
        WG->node[i].kout=0;
        WG->node[i].kout=0;
        WG->node[i].sin=0;
        WG->node[i].sout=0;
        WG->node[i].mem_in=1;
        WG->node[i].mem_out=1;
        WG->node[i].out=calloc(1,sizeof(int));
        WG->node[i].out[0]=-1;
        WG->node[i].in=calloc(1,sizeof(int));
        WG->node[i].in[0]=-1;
        WG->node[i].w_in=calloc(1,sizeof(int));
        WG->node[i].w_out=calloc(1,sizeof(int));
    }
    return WG;
}

void w_graph_free_nodes(W_GRAPH* WG, int N_nodes){
    int i;
    for(i=0;i<N_nodes;i++)
    {
        free(WG->node[i].out);
        free(WG->node[i].in);
		free(WG->node[i].w_out);
		free(WG->node[i].w_in);
    }
    free(WG->node);
    return;
}

/****************************************************************************
 * Add edges *
 ****************************************************************************/
/// warning: Self-loops accepted ///
void w_graph_add_multi_link(W_GRAPH * WG, int N_nodes, int origin, int dest, int weight){
	int neigh,dummy;
	if(weight > 0)
	{
		WG->node[origin].sout +=weight;
		WG->node[dest].sin	 +=weight;
		neigh = find_value_int(WG->node[origin].out, dest, WG->node[origin].kout);
		if(neigh<0)// new connection
		{
			WG->node[origin].kout++;
			if(WG->node[origin].kout > WG->node[origin].mem_out)
			{
				dummy = WG->node[origin].mem_out;
				WG->node[origin].mem_out = dummy*2;
				WG->node[origin].out		= safe_int_realloc(WG->node[origin].out, dummy, WG->node[origin].mem_out, -1);
				//safe_int_realloc(WG->node[origin].out,dummy,node[origin].mem_out);
				WG->node[origin].w_out	= safe_int_realloc(WG->node[origin].w_out, dummy, WG->node[origin].mem_out, -1);
				//safe_int_realloc(WG->node[origin].w_out,dummy,node[origin].mem_out);
			}
			WG->node[origin].out[WG->node[origin].kout-1] 	= dest;
			WG->node[origin].w_out[WG->node[origin].kout-1] 	= weight;
		}else{ // existing connection
			WG->node[origin].w_out[neigh]					+=	weight;
		}
		neigh = find_value_int(WG->node[dest].in, origin, WG->node[dest].kin);
		if(neigh<0)// new connection
		{
			WG->node[dest].kin++;
			if(WG->node[dest].kin > WG->node[dest].mem_in)
			{
				dummy 					= WG->node[dest].mem_in;
				WG->node[dest].mem_in 	= 2*dummy;
				WG->node[dest].in 		= safe_int_realloc(WG->node[dest].in, dummy, WG->node[dest].mem_in, -1);
				//safe_int_realloc(WG->node[dest].in,dummy,node[dest].mem_in);
				WG->node[dest].w_in 		= safe_int_realloc(WG->node[dest].w_in, dummy, WG->node[dest].mem_in, -1);
				//safe_int_realloc(WG->node[dest].w_in,dummy,node[dest].mem_in);
			}
			WG->node[dest].in[WG->node[dest].kin-1]	= origin;
			WG->node[dest].w_in[WG->node[dest].kin-1]	= weight;
		}else{ // existing connection
			WG->node[dest].w_in[neigh] 				+= weight;
		}
	}
    return;
}


void w_graph_add_multi_link_undirected(W_GRAPH * WG, int N_nodes, int origin, int dest, int weight){
    if(weight > 0)
    {
		// keep sin-sout for compatibility issues and add twice if self_loop
		// update s
		WG->node[origin].sout	+= weight;
		WG->node[origin].sin		+= weight;
		if(dest!=origin)
		{
			WG->node[dest].sout 	+= weight;
			WG->node[dest].sin 	+= weight;
		}
		// check new connection
		int neigh,dummy;
		neigh = find_value_int(WG->node[origin].out, dest, WG->node[origin].kout);
		if(neigh<0)// new connection
		{
		// update k
			WG->node[origin].kout++;
			WG->node[origin].kin++;
		// if not enough space, allocate more
			if(WG->node[origin].kout>WG->node[origin].mem_out)
			{
				dummy 					= WG->node[origin].mem_out;
				WG->node[origin].mem_out	= 2*dummy;
				WG->node[origin].out		= safe_int_realloc(WG->node[origin].out, dummy, WG->node[origin].mem_out, -1);
				//safe_int_realloc(WG->node[origin].out,dummy,node[origin].mem_out);
				//printf("lalal am here , node: %d, dest: %d, k:%d, mem:%d\n", origin, dest, node[origin].kout, node[origin].mem_out); 				fflush(stdout);
				//safe_int_realloc(WG->node[origin].w_out,dummy,node[origin].mem_out);
				WG->node[origin].w_out	= safe_int_realloc(WG->node[origin].w_out, dummy, WG->node[origin].mem_out, -1);
			}
			// store neighbor
			WG->node[origin].out[WG->node[origin].kout-1] = dest;
			WG->node[origin].w_out[WG->node[origin].kout-1]		= weight;
			if(dest!=origin)
			{
				WG->node[dest].kout++;
				WG->node[dest].kin++;
				if(WG->node[dest].kout > WG->node[dest].mem_out)
				{
					dummy 					= WG->node[dest].mem_out;
					WG->node[dest].mem_out	= 2*dummy;
					WG->node[dest].out		= safe_int_realloc(WG->node[dest].out, dummy, WG->node[dest].mem_out,-1);
					//safe_int_realloc(WG->node[dest].out,dummy,node[dest].mem_out);
					WG->node[dest].w_out		= safe_int_realloc(WG->node[dest].w_out, dummy, WG->node[dest].mem_out,-1);
					//safe_int_realloc(WG->node[dest].w_out,dummy,node[dest].mem_out);
				}
				WG->node[dest].out[WG->node[dest].kout-1]   = origin;
				WG->node[dest].w_out[WG->node[dest].kout-1] = weight;
			}
		}else{ // existing connection
			WG->node[origin].w_out[neigh]   += weight;
			if(origin != dest)
			{
				neigh = find_value_int(WG->node[dest].out, origin, WG->node[dest].kout);
				assert(neigh>=0);
				WG->node[dest].w_out[neigh] += weight;
			}
		}
	}
    return;
}

/****************************************************************************
 * Printing *
 ****************************************************************************/
void w_graph_print_adj_list(W_GRAPH* WG, int N_nodes, char* output){
    int i,j;
    FILE* fil=open_file("w", output);
    fprintf(fil,"## Adjancenncy list (directed weighted format): Node id_source Node_id target weight \n");fflush(stdout);
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<WG->node[i].kout;j++)
        {
            fprintf(fil,"%d %d %d\n",WG->node[i].idnum, WG->node[WG->node[i].out[j]].idnum, WG->node[i].w_out[j]);
        }
    }
    fclose(fil);    
    return;
}
/****************************************************************************
 * S and K's *
 ****************************************************************************/

int ** w_graph_compute_s(W_GRAPH* WG, int N_nodes){
    int** s=cast_mat_int(2,N_nodes);
    int i;
    for(i=0;i<N_nodes;i++)
    {
        s[0][i] = WG->node[i].sout;
        s[1][i] = WG->node[i].sin;
    }
    return s;
}


int ** w_graph_compute_k(W_GRAPH* WG, int N_nodes){
    int** k=cast_mat_int(2,N_nodes);
    int i;
    for(i=0;i<N_nodes;i++)
    {
        k[0][i]=WG->node[i].kout;
        k[1][i]=WG->node[i].kin;
    }
    return k;
}


double ** w_graph_compute_k_analitic(W_GRAPH* WG, int N_nodes, int self_opt){
    double** k=cast_mat_double(2,N_nodes);
    int **s=w_graph_compute_s(WG, N_nodes);
    long int T=sum_vec_int(s[0],N_nodes);
    int i,j;
    for(i=0;i<N_nodes;i++)
    {
        k[0][i]=(double)N_nodes;
        k[1][i]=(double)N_nodes;
        if(self_opt<=0)
        {
            k[0][i]-=1;
            k[1][i]-=1;

        }
        for(j=0;j<N_nodes;j++)
        {
            if(j!=i)
            {
                k[0][i]-=exp(-(double)WG->node[i].sout*((double)WG->node[j].sin/(long double)T));
                k[1][i]-=exp(-(double)WG->node[i].sin*((double)WG->node[j].sout/(long double)T));
            }else{
            	if(self_opt>0)
				{
	                k[0][i]-=exp(-(double)WG->node[i].sout*((double)WG->node[j].sin/(long double)T));
	                k[1][i]-=exp(-(double)WG->node[i].sin*((double)WG->node[i].sout/(long double)T));					
				}
            }
        }
    }
    free_mat_int(s,2);
    return k;
}



double ** w_graph_compute_k_analitic_from_s_directed(int** s, int N_nodes, int self_opt){
    double** k=cast_mat_double(4,N_nodes); // 2 firts out-in, then sigma.
    long int T=(long int)sum_vec_int(s[0],N_nodes);
    int i,j;
    for(i=0;i<N_nodes;i++)
    {
        k[0][i]=(double)N_nodes;
        k[1][i]=(double)N_nodes;
        k[2][i]=(double)N_nodes;
        k[3][i]=(double)N_nodes;
        if(self_opt<=0)
        {
            k[0][i]-=1;
            k[1][i]-=1;
            k[2][i]-=1;
            k[3][i]-=1;
        }
        for(j=0;j<N_nodes;j++)
        {
            if(j!=i)
            {
                k[0][i]-=exp(-(double)s[0][i]*((double)s[1][j]/(long double)T));
                k[1][i]-=exp(-(double)s[1][i]*((double)s[0][j]/(long double)T));
				k[2][i]-=exp(-2*(double)s[0][i]*((double)s[1][j]/(long double)T));
				k[3][i]-=exp(-2*(double)s[1][i]*((double)s[0][j]/(long double)T));
            }else{
				if(self_opt>0)
				{
		    		k[0][i]-=exp(-(double)s[0][i]*((double)s[1][i]/(long double)T));
		    		k[1][i]-=exp(-(double)s[1][i]*((double)s[0][i]/(long double)T));
		    		k[2][i]-=exp(-2*(double)s[0][i]*((double)s[1][i]/(long double)T));
		    		k[3][i]-=exp(-2*(double)s[1][i]*((double)s[0][i]/(long double)T));

				}
	    	}
        }
		k[2][i] = k[2][i]-k[0][i]; // sigma = N- k - exp 2!
		k[3][i] = k[3][i]-k[1][i];
    }
    return k;
}

double ** w_graph_compute_k_analitic_from_s_undirected(int* s, int N_nodes, int self_opt){
    double** k=cast_mat_double(2,N_nodes);
    long int T=(long int)sum_vec_int(s,N_nodes);
    int i,j;
    for(i=0;i<N_nodes;i++)
    {
        k[0][i]=(double)N_nodes;
		k[1][i]=(double)N_nodes;
        if(self_opt<=0)
		{
			k[0][i]-=1;		
			k[1][i]-=1;		
		} 
        for(j=0;j<N_nodes;j++)
        {
            if(j!=i)
            {
                k[0][i]-=exp(-(double)s[i]*((double)s[j]/(long double)T));
                k[1][i]-=exp(-2*(double)s[i]*((double)s[j]/(long double)T));

            }else{
				if(self_opt>0)
				{
		    		k[0][i]-=exp(-(double)s[i]*((double)s[j]/(long double)T));
		    		k[1][i]-=exp(-2*(double)s[i]*((double)s[j]/(long double)T));
				}
	    	}
        }
		k[1][i] = k[1][i]-k[0][i];
    }
    return k;
}



/****************************************************************************
 * Snn and Knn's *
 ****************************************************************************/
double ** w_graph_compute_s_nn(W_GRAPH* WG, int N_nodes, int weight, int opt_dir){
    double ** s_nn=cast_mat_double(4, N_nodes);
    int i,j;
    //FILE* prova=open_file("a","check.dat");
    for(i=0;i<N_nodes;i++)
    {
        s_nn[0][i]=s_nn[1][i]=s_nn[2][i]=s_nn[3][i]=0;
        if(WG->node[i].kout>0)
	{
	for(j=0;j<WG->node[i].kout;j++)
        {
            if(weight>0)
            {
                s_nn[1][i]+=(double)(WG->node[i].w_out[j])*(double)WG->node[WG->node[i].out[j]].sin; // average out-neigh, in strength
                s_nn[0][i]+=(double)(WG->node[i].w_out[j])*(double)WG->node[WG->node[i].out[j]].sout; // average out-neigh, out
            }else if(weight==0){
                s_nn[1][i]+=(double)(WG->node[i].w_out[j])*(double)WG->node[WG->node[i].out[j]].kin; // weighted degree
                s_nn[0][i]+=(double)(WG->node[i].w_out[j])*(double)WG->node[WG->node[i].out[j]].kout; // weighted degree
            }else{
                s_nn[1][i]+=(double)WG->node[WG->node[i].out[j]].kin; // average out-neigh, in strength
                s_nn[0][i]+=(double)WG->node[WG->node[i].out[j]].kout; // average out-neigh, out                                
            }
        }
    	}
	if((WG->node[i].kin>0) && (opt_dir>0))
	{
        for(j=0;j<WG->node[i].kin;j++)
        {        
            if(weight>0)
            {
                s_nn[3][i]+=WG->node[i].w_in[j]*(double)WG->node[WG->node[i].in[j]].sin; // average in-neigh, in
                s_nn[2][i]+=WG->node[i].w_in[j]*(double)WG->node[WG->node[i].in[j]].sout; // average in-neigh, out
            }else if(weight==0){
                s_nn[3][i]+=WG->node[i].w_in[j]*(double)WG->node[WG->node[i].in[j]].kin; // average in-neigh, in
                s_nn[2][i]+=WG->node[i].w_in[j]*(double)WG->node[WG->node[i].in[j]].kout; // average in-neigh, out
            }else{
                s_nn[3][i]+=(double)WG->node[WG->node[i].in[j]].kin; // average out-neigh, in strength
                s_nn[2][i]+=(double)WG->node[WG->node[i].in[j]].kout; // average out-neigh, out                
            }

        }
	}
        if(weight>=0)
        {
	    if(WG->node[i].kout>0)
	    {
            s_nn[0][i]/=(double)WG->node[i].sout;
            s_nn[1][i]/=(double)WG->node[i].sout;
	    }
	    if((WG->node[i].kin>0) && (opt_dir>0))
	    {
	    	s_nn[2][i]/=(double)WG->node[i].sin;
            s_nn[3][i]/=(double)WG->node[i].sin;
	    }
        //fprintf(prova,"%d %d %lf\n",i,node[i].sout,s_nn[0][i]);
        }else{
	    if(WG->node[i].kout>0)
	    {
            s_nn[0][i]/=(double)WG->node[i].kout;
            s_nn[1][i]/=(double)WG->node[i].kout;
	    }
	    if((WG->node[i].kin>0) && (opt_dir>0))
	    {
            s_nn[2][i]/=(double)WG->node[i].kin;
            s_nn[3][i]/=(double)WG->node[i].kin;                
	    }
        }
    }
    //fclose(prova);
    return s_nn;
}

/****************************************************************************
 * Clustering *
 ****************************************************************************/

double ** w_graph_compute_clust(W_GRAPH * WG, int N_nodes){ // 2 cols: unweighted, weighted
    //printf("#### Computing c ####\n"); fflush(stdout);
    int i,j,k,l;
    int c_glo, c_glo_den,e;
    double ee;
    c_glo=c_glo_den=0;
    double** c=cast_mat_double(2,N_nodes);
    ///for per cada node i
    for(i=0;i<N_nodes;++i)
    {
	if(WG->node[i].kout>0)
	{
    //printf("start, %d, %d \n",i,node[i].kin);fflush(stdout);
    	e=ee=0;
	//printf("node=%i\n",i);
	///for per cada vei j
	for(j=0;j<WG->node[i].kout;++j)
	{
	    //printf("vei=%i amb degree=%i\n",node[i].out[j],node[WG->node[i].out[j]].k);
	    ///for per els altres nodes veins de i. comencem a j+1 per no repetir
	    for(k=j+1;k<WG->node[i].kout;++k)
	    {
		//printf("segon_vei=%i\n",node[i].out[k]);
		///un for per veure els veins del node j i veure si esta conectat amb el node k
		for (l=0;l<WG->node[WG->node[i].out[j]].kout;++l)
		{
		    //printf("l=%i\n",l);
		    if(WG->node[WG->node[i].out[j]].out[l]==WG->node[i].out[k]) 
		    {
			e=e+1;
			ee+=(WG->node[i].w_out[j]+WG->node[i].w_out[k]);
		    } 
		    //printf("   node=%i vei=%i veivei=%i veicomparar=%i  e=%i\n",i,node[i].out[j],node[WG->node[i].out[j]].out[l],node[i].out[k],e);
		}
	    }
	}                        
	if(e!=0)
        {
	    c_glo=c_glo+e; // triangles
	    c_glo_den=c_glo_den+WG->node[i].kout*(WG->node[i].kout-1)/2; // pairs of neighbours
	    //cc=cc+(double)e*2./node[i].k/(WG->node[i].k-1); // l
	    c[0][i]= (double)e*2./WG->node[i].kout/(WG->node[i].kout-1); // clustering
	    c[1][i]= (double)ee/(2.*WG->node[i].sout)/(WG->node[i].kout-1); // weighted clust
	}
    	}
        //printf("stop, %d\n",i);fflush(stdout);
    }
    ///and we take the mean
    //cc=cc/(n-number_k[1]);
    //printf("global clustering= %f\n",(double)c_glo/(double)c_glo_den);
    return c;
}



/****************************************************************************
 * Weight funcs *
 ****************************************************************************/
int w_graph_total_weight( W_GRAPH* WG, int N_nodes){
    int i;
    int T;
    T=0;
    for(i=0;i<N_nodes;i++)
    {
        T+=(int)WG->node[i].sout;
    }
    return T;
}

int w_graph_total_edges( W_GRAPH* WG, int N_nodes){
    // counts self loops properly (bit more slower... but ok)
    int i,j;
    int T,aux2;
    T=aux2=0;
    for(i=0;i<N_nodes;i++)
    {
        T+=WG->node[i].kout;
        for(j=0;j<WG->node[i].kout;j++)
        {
            aux2++;
            if(WG->node[i].out[j]==i)
            {
                aux2++;
            }
        }
    }
    //return T;
    return aux2;
}


int * w_graph_compute_w(W_GRAPH* WG, int N_nodes, int* aux, int zeros){
    //zeros not implemented
    int E = w_graph_total_edges( WG, N_nodes);
    int* w=cast_vec_int(E);
    int i,j,aux2,mem;
    mem=E;
    aux2=0;
    if(zeros>0)
    {
        for(i=0;i<N_nodes;i++)
        {
            for(j=0;j<WG->node[i].kout;j++)
            {
                if(aux2+1>mem)
                {
                    w=realloc(w, sizeof(int)*2*mem);
                    mem=2*mem;
                }
                w[aux2]=(int)WG->node[i].w_out[j];
                aux2++;
                if(WG->node[i].out[j]==i) // count twice
                {
                    if(aux2+1>mem)
                    {
                        w=realloc(w, sizeof(int)*2*mem);
                        mem=2*mem;
                    }
                    w[aux2]=(int)WG->node[i].w_out[j];
                    aux2++;
                }

            }
        }
    }else{
        for(i=0;i<N_nodes;i++)
        {
            for(j=0;j<WG->node[i].kout;j++)
            {
                if(aux2+1>mem)
                {
                    w=realloc(w, sizeof(int)*2*mem);
                    mem=2*mem;
                }
                w[aux2]=(int)WG->node[i].w_out[j];
                aux2++;
                if(WG->node[i].out[j]==i)
                {
                    if(aux2+1>mem)
                    {
                        w=realloc(w, sizeof(int)*2*mem);
                        mem=2*mem;
                    }
                    w[aux2]=(int)WG->node[i].w_out[j];
                    aux2++;
                }
                
            }
        }
    }
    //printf("E: %d, aux2:%d delta:%d", E, aux2, E-aux2); fflush(stdout);
    *aux=aux2;
    realloc(w,sizeof(int)*aux2); // fix final size
    return w;
}

double** w_graph_compute_p_w_analitic_from_s_undirected(int maxt, double binn, int* s, int N_nodes, int self_opt, int* len){
    // computes already normalized p
    printf("Computing w, tmax:%d\n",maxt); fflush(stdout);
    int i,j,t,aux;
    int T=sum_vec_int(s,N_nodes);
    double norm,mu;
    double** p=cast_mat_double(2,N_nodes);
    mu=0;
    norm=0;
    aux=0;
    t=1;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<i;j++)
        {
            norm+=1.-(double)exp(-(double)s[i]*(double)s[j]/(double)T);
        }
        if(self_opt>0)
        {
            norm+=1.-(double)exp(-(double)s[i]*(double)s[i]/(double)T);
        }
    }
    while(t<maxt+1)
    {
        p[0][aux]=t;
        p[1][aux]=0;
        //fact=factorial(t);
        for(i=0;i<N_nodes;i++)
        {
            for(j=0;j<i;j++)
            {
                mu = (double)s[i]*(double)s[j]/(double)T;
                p[1][aux]+=gsl_ran_poisson_pdf (t, mu);
                if(t==0)norm+=(double)exp(-mu);
            }
            if(self_opt>0)
            {
                mu = (double)s[i]*(double)s[i]/(double)T;
                if(t==0)norm+=(double)exp(-mu);
                p[1][aux]+=gsl_ran_poisson_pdf (t, mu);
            }
        }
        p[1][aux]=p[1][aux]/norm;
        //printf("done %d\n",t); fflush(stdout);
        if(t<10){
            t++;
        }else{
            t=t*binn;
        }
        aux++;
    }
    (*len)=aux;
    printf("done: %d bins\n",aux); fflush(stdout);
    realloc(p[0],sizeof(double)*aux);
    realloc(p[1],sizeof(double)*aux);
    return p;
}

double** w_graph_compute_p_w_analitic_from_s_directed(int maxt, double binn, int** s, int N_nodes, int self_opt, int* len){
    // computes already normalized p
    //printf("Computeing w, tmax:%d\n",maxt); fflush(stdout);
    printf("Computing w, tmax:%d\n",maxt); fflush(stdout);
    int i,j,t,aux;
    int T=sum_vec_int(s[0],N_nodes);
    double norm,mu;
    double** p=cast_mat_double(2,N_nodes);
    mu=0;
    norm=(double)N_nodes;
    aux=0;
    t=1;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<N_nodes;j++)
        {
            mu = (double)s[0][i]*(double)s[1][j]/(double)T;
            if(i!=j)norm+=1.-(double)exp(-mu);
        }
        if(self_opt>0)
        {
            mu = (double)s[0][i]*(double)s[1][i]/(double)T;
           norm+=1.-(double)exp(-mu);
        }
    }
    while(t<maxt+1)
    {
        p[0][aux]=t;
        p[1][aux]=0;
        for(i=0;i<N_nodes;i++)
        {
            for(j=0;j<N_nodes;j++)
            {
                if(j!=i)
                {
                    mu = (double)s[0][i]*(double)s[1][j]/(double)T;
                    p[1][aux]+=gsl_ran_poisson_pdf (t, mu);
                }
            }
            if(self_opt>0)
            {
                mu = (double)s[0][i]*(double)s[1][j]/(double)T;
                p[1][aux]+=gsl_ran_poisson_pdf (t, mu);
            }
        }
        p[1][aux]=p[1][aux]/norm;
        //printf("done %d\n",t); fflush(stdout);
        if(t<10){
            t++;
        }else{
            t=t*binn;
        }
        aux++;
    }
    (*len)=aux;
    //printf("done\n"); fflush(stdout);
    realloc(p[0],sizeof(double)*aux);
    realloc(p[1],sizeof(double)*aux);
    return p;
}

double * w_graph_compute_w_ss(W_GRAPH* WG, int N_nodes, int weight){
//int * w_graph_compute_w_ss(W_GRAPH* node, int N_nodes, int weight){
    int E;
    E=w_graph_total_edges(WG,N_nodes);
    //int* ww=cast_vec_int(E);
    double* ww=cast_vec_double(E);
    int i,j,aux;
    aux=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<WG->node[i].kout;j++)
        {
            if (weight>0)
            {
                ww[aux]=((double)WG->node[i].sout)*((double)WG->node[WG->node[i].out[j]].sin);
            }else{
                ww[aux]=((double)WG->node[i].kout)*((double)WG->node[WG->node[i].out[j]].kin);
            }
            aux++;
            if(WG->node[i].out[j]==i) // if self loop, count twice
            {
                if (weight>0)
                {
                    ww[aux]=((double)WG->node[i].sout)*((double)WG->node[WG->node[i].out[j]].sin);
                }else{
                    ww[aux]=((double)WG->node[i].kout)*((double)WG->node[WG->node[i].out[j]].kin);
                }
                aux++;
            }
        }
    }
    assert(aux==E);
    return ww;
}


//double ** w_graph_compute_xy(w_graph * node, int N_nodes){
    
    //int i,j;
    //double** y2=cast_mat_double(3,N_nodes);
    //for(i=0;i<N_nodes;i++)
    //{
	//if(WG->node[i].kout>0)
	//{
        //y2[0][i]=0;
	    //for(j=0;j<node[i].kout;j++)
	    //{
	    	//y2[0][i]+=((double)node[i].w_out[j])*((double)node[i].w_out[j]); // x
	    //}
	//y2[1][i]=node[i].sout*node[i].sout; // y
	//y2[2][i]=y2[1][i]*y2[0][i]; // xy
    	//}
    //}
    //return y2;
//}


double ** w_graph_compute_Y2(W_GRAPH * WG, int N_nodes, int opt_dir){
    
    int i,j;
    double** y2=cast_mat_double(2,N_nodes);
    for(i=0;i<N_nodes;i++)
    {
	y2[0][i]=0;
	if(WG->node[i].kout>0)
	{
	    for(j=0;j<WG->node[i].kout;j++)
	    {
	    	y2[0][i]+=((double)WG->node[i].w_out[j])*((double)WG->node[i].w_out[j]);
/*            if((WG->node[i].w_out[j]<1) )
            {
                printf("%f\n",(double)WG->node[i].w_out[j]);
            }
 */
	    }
	    y2[0][i]/=(((double)WG->node[i].sout)*((double)WG->node[i].sout));
    }
	y2[1][i]=0;
	if((WG->node[i].kin>0) && (opt_dir>0))
	{        
	    for(j=0;j<WG->node[i].kin;j++)
	    {
            	y2[1][i]+=((double)WG->node[i].w_in[j])*((double)WG->node[i].w_in[j]);            
	    }
	    y2[1][i]/=(((double)WG->node[i].sin)*((double)WG->node[i].sin));
    }
/*    if((y2[0][i] >0)|( y2[1][i] >0))
       {
           printf("out y2:%f \t s: %d \t k: %d\n in y2:%f \t s: %d \t k: %d\n\n",y2[0][i],WG.WG->node[i].sout,WG->node[i].kout,y2[1][i],WG->node[i].sin,WG->node[i].kin); fflush(stdout);
       }
*/
    }
    return y2;
}

/****************************************************************************
 * Entropy & LIkelyhood funcs *
 ****************************************************************************/
double w_graph_entropy(W_GRAPH* WG, int N_nodes){
	double S,p;
	int i,j,t;
	int T = w_graph_total_weight(WG,N_nodes);
	S = 0;
	for(i=0;i<N_nodes;i++) // all edges
	{
		for(j=0;j<WG->node[i].kout;j++)
		{
			t = WG->node[i].w_out[j];
			p = (double)t/(double)T;
			S+= p*log(p);
		}		
	}
	return -S;
}


void w_graph_print_entropy(double* seq,int len,char* output){
	int bins;
	double mins,maxs,eps;
	bins = (int)len/20.;
	if(bins<5)
	{
		bins=5;
	}else if(bins>50){
		bins=50;
	}
	mins = min_value_double(seq,len);
	maxs = max_value_double(seq,len);
	eps = (maxs-mins)/100.;
    gsl_histogram* h1=histogram_double(seq,mins-mins*eps,maxs+maxs*eps,bins,len);
	print_acc(output,h1,h1);
	return;
}

double w_graph_loglikelyhood(W_GRAPH* WG,int N_nodes,double** wij){
	double L=0;
	int i,j,t;
	double p,mu;
	//int T = w_graph_total_weight(WG->node,N_nodes);
	for(i=0;i<N_nodes;i++) // all edges
	{
		for(j=0;j<WG->node[i].kout;j++)
		{
			t = WG->node[i].w_out[j];
			mu = wij[i][WG->node[i].out[j]];
			if(mu>0)
			{
				p = gsl_ran_poisson_pdf (t, mu);
				//printf("Mu:%f p:%f t:%d",mu,p,t);fflush(stdout);
				L+= log(p);
			}
		}
				
	}
	return L;
}
/****************************************************************************
 * aLL STATS *
 ****************************************************************************/
void w_graph_node_stats_list(W_GRAPH* WG, int N_nodes, int run, double av_k, int opt_dir, int opt_clust, int self_opt){
    //p(k),p(s),p(w),w(sin sout), s_nn, k_nn, s(k)
    int **k=w_graph_compute_k(WG, N_nodes);
    int **s=w_graph_compute_s(WG, N_nodes);
    //int E;
    //long int T=sum_vec_int(s[0],N_nodes);
    //int *wkk2=w_graph_compute_w_ss(WG->node, N_nodes, 0);
    double ** ss_n=w_graph_compute_s_nn(WG, N_nodes,1, opt_dir);// 4 rows
    double ** kk_n=w_graph_compute_s_nn(WG, N_nodes,-1, opt_dir);// 4 rows
    double ** kkw_n=w_graph_compute_s_nn(WG, N_nodes,0, opt_dir);// 4 rows
    double **y2 = w_graph_compute_Y2(WG, N_nodes, opt_dir); // 2 rows
    char cadena[100];
    int i; 
    // Get analytical predictions //
    // k(s) //
	if(opt_dir>0)
	{
		sprintf(cadena,"N%davs%8.5fnode_list.list",N_nodes,av_k);
	}else{
		sprintf(cadena,"N%davs%8.5f_undir_node_list.list",N_nodes,av_k);
	}
    
    FILE* fil=open_file("w", cadena);
    if(opt_dir==1)
    {
    	double ** k_anal=w_graph_compute_k_analitic(WG, N_nodes, self_opt);
		// Note for indices on s_nn:
		// [0] -> \sum t_{ij} s_out_j / s_out_i
		// [1] -> \sum t_{ij} s_in_j / s_out_i
		// [2] -> \sum t_{ji} s_out_j / s_in_i
		// [3] -> \sum t_{ji} s_in_j / s_in_i
    	fprintf(fil,"# Node_num\tk\tk_anal\ts\tY2\tk_nn\tk^w_nn\ts^w_nn (out) then (in) # \n");
    	for(i=0;i<N_nodes;i++)
    	{
        	fprintf(fil,"%d %d %.3f %d %f %f %f %f %d %.3f %d %f %f %f %f\n",i,k[0][i],k_anal[0][i],s[0][i],
			y2[0][i],kk_n[1][i],kkw_n[1][i],ss_n[1][i],k[1][i],k_anal[1][i],s[1][i],y2[1][i],kk_n[2][i],kkw_n[2][i],ss_n[2][i]);
    	}
	free_mat_double(k_anal,2);
    }else{
	double ** k_anal=w_graph_compute_k_analitic(WG, N_nodes, self_opt);
	if(opt_clust==1)
	{
	    double ** c= w_graph_compute_clust(WG, N_nodes); // 2 rows
	    fprintf(fil,"# Node_num\tk\tk_anal\ts\tY2\tk_nn\tk^w_nn\ts^w_nn\tc\tc^w# \n");
	    for(i=0;i<N_nodes;i++)
	    {
        	fprintf(fil,"%d %d %.3f %d %f %f %f %f %f %f\n",i,k[0][i],k_anal[0][i],s[0][i],y2[0][i],kk_n[0][i],
			kkw_n[0][i],ss_n[0][i],c[0][i],c[1][i]);
	    }
	    free_mat_double(c,2);
	}else{
	    fprintf(fil,"# Node_num\tk\tk_anal\ts\tY2\tk_nn\tk^w_nn\ts^w_nn\n");
	    for(i=0;i<N_nodes;i++)
	    {
        	fprintf(fil,"%d %d %.3f %d %f %f %f %f\n",i,k[0][i],k_anal[0][i],s[0][i],y2[0][i],kk_n[0][i],
			kkw_n[0][i],ss_n[0][i]);
	    }
	free_mat_double(k_anal,2);
	}
    }
    fclose(fil);
    free_mat_int(s,2);
    free_mat_int(k,2);
    free_mat_double(ss_n,4);
    free_mat_double(kk_n,4);
    free_mat_double(kkw_n,4);
    free_mat_double(y2,2);
    return;
}


void w_graph_all_stats(W_GRAPH* WG, int N_nodes, int run, double bin_exp, double av_k, int opt_dir, int self_opt, int w_anal){
    //w(sin sout), w(kin,kout), w
    int **k=w_graph_compute_k(WG, N_nodes);
    int **s=w_graph_compute_s(WG, N_nodes);
    int E;
    int *w=w_graph_compute_w(WG, N_nodes, &E, -1);
    double *wss=w_graph_compute_w_ss(WG, N_nodes, 1);
    double *wkk=w_graph_compute_w_ss(WG, N_nodes, -1);
    char cadena[100];
    
    gsl_histogram* h1;
    
    double* sout;
    double* xranges;
    int xbins;
    double** yy;
    
    sout=vec_int_to_double(w,E);
    int q=max_value_int(w,E);
    h1=histogram_double(sout,0,q,q,E);
    free(sout);
    //sprintf(cadena,"run_%dN%d_w.hist",run,N_nodes);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%davs%8.5f_w.hist",N_nodes,av_k);
	}else{
	    sprintf(cadena,"N%davs%8.5f_undir_w.hist",N_nodes,av_k);
	}

    print_acc(cadena, h1, h1);
    gsl_histogram_free(h1);

    /// analitical w /////
    if(w_anal>0)
    {
        double** pp;
        int lenn;
        if(opt_dir>0)
        {
            pp = w_graph_compute_p_w_analitic_from_s_directed(10*q,1.5,s,N_nodes, self_opt, &lenn);
			sprintf(cadena,"N%davs%8.5f_w_anal.hist",N_nodes,av_k);
        }else{
            pp = w_graph_compute_p_w_analitic_from_s_undirected(10*q,1.5,s[0],N_nodes, self_opt, &lenn);
			sprintf(cadena,"N%davs%8.5f_undir_w_anal.hist",N_nodes,av_k);
        }
        FILE* fil=open_file("w", cadena);
        fprintf(fil,"# t p(t) # \n");
        int i;
        for(i=0;i<lenn;i++)
        {
            fprintf(fil,"%.8f %.8f\n",pp[0][i],pp[1][i]);
        }
        fclose(fil);
        //free(pp[0]);
        //free(pp[1]);
        //free(pp);
    }
    sout=vec_int_to_double(w,E);
    xranges=log_bins_double(0, max_value_double(wss,E) , bin_exp, &xbins);
    yy=y_of_x(wss, sout, xranges,  E,  xbins);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%davs%8.5f_w_s_oi.hist",N_nodes,av_k);
	}else{
	    sprintf(cadena,"N%davs%8.5f_undir_w_s_oi.hist",N_nodes,av_k);
	}
    print_hist2d_mean(cadena, yy[1], yy[2], yy[0], xbins-1);
    free(sout);
    free(wss);
    free(xranges);
    free_mat_double(yy,4);

    sout=vec_int_to_double(w,E);
    xranges=log_bins_double(0, max_value_double(wkk,E) , bin_exp, &xbins);
    yy=y_of_x(wkk, sout, xranges,  E,  xbins);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%davs%8.5f_w_k_oi.hist",N_nodes,av_k);
	}else{
	    sprintf(cadena,"N%davs%8.5f_undir_w_k_oi.hist",N_nodes,av_k);
	}

    print_hist2d_mean(cadena, yy[1], yy[2], yy[0], xbins-1);
    free(sout);
    free(wkk);
    free(xranges);
    free_mat_double(yy,4);

    
//// free all
    free_mat_int(s,2);
    free_mat_int(k,2);
    free(w);
    return;
}


/****************************************************************************
 * Ensembles averaging stats *
 ****************************************************************************/

// Nodes
void w_graph_node_stats_ensemble(W_GRAPH* WG, int N_nodes, double** container, double ** container2, int** node_nonzero ,double* T_container, int opt_dir, int opt_clust ){
    //p(k),p(s),p(w),w(sin sout), s_nn, k_nn, s(k)
    int **k=w_graph_compute_k(WG, N_nodes);
    int **s=w_graph_compute_s(WG, N_nodes);
    //int E;
    int T=sum_vec_int(s[0],N_nodes);
    //int *wkk2=w_graph_compute_w_ss(WG->node, N_nodes, 0);
    double ** ss_n=w_graph_compute_s_nn(WG, N_nodes,1, opt_dir);// 4 rows
    double ** kk_n=w_graph_compute_s_nn(WG, N_nodes,-1, opt_dir);// 4 rows
    double ** kkw_n=w_graph_compute_s_nn(WG, N_nodes,0, opt_dir);// 4 rows
    double **y2 = w_graph_compute_Y2(WG, N_nodes, opt_dir); // 2 rows
    //double **xy = w_graph_compute_xy(WG->node, N_nodes); // 2 rows
    //char cadena[100];
		// Note for indices on s_nn:
		// [0] -> \sum t_{ij} s_out_j / s_out_i
		// [1] -> \sum t_{ij} s_in_j / s_out_i
		// [2] -> \sum t_{ji} s_out_j / s_in_i
		// [3] -> \sum t_{ji} s_in_j / s_in_i

    int i;
    if(opt_dir==1)
    { 
    	for(i=0;i<N_nodes;i++)
    	{
	    if(s[0][i] > 0)
	    {
		//assert(s[0][i]>=k[0][i]);
		node_nonzero[i][0]+=1;
		container[i][0]+=(double)k[0][i];
		container2[i][0]+=(double)k[0][i]*k[0][i];
		container[i][2]+=(double)s[0][i];
		container2[i][2]+=(double)s[0][i]*s[0][i];
		container[i][3]+=(double)y2[0][i];
		container2[i][3]+=(double)y2[0][i]*y2[0][i];
		container[i][4]+=(double)kk_n[1][i];
		container2[i][4]+=(double)kk_n[1][i]*kk_n[1][i];
		container[i][5]+=(double)kkw_n[1][i];
		container2[i][5]+=(double)kkw_n[1][i]*kkw_n[1][i];
		container[i][6]+=(double)ss_n[1][i];
		container2[i][6]+=(double)ss_n[1][i]*ss_n[1][i];

	    }
	    if(s[1][i] > 0)
	    {
		//assert(s[1][i]>=k[1][i]);
		node_nonzero[i][1]+=1;
		container[i][7]+=(double)k[1][i];
		container2[i][7]+=(double)k[1][i]*k[1][i];
		container[i][9]+=(double)s[1][i];
		container2[i][9]+=(double)s[1][i]*s[1][i];
		container[i][10]+=(double)y2[1][i];
		container2[i][10]+=(double)y2[1][i]*y2[1][i];
		container[i][11]+=(double)kk_n[2][i];
		container2[i][11]+=(double)kk_n[2][i]*kk_n[2][i];
		container[i][12]+=(double)kkw_n[2][i];
		container2[i][12]+=(double)kkw_n[2][i]*kkw_n[2][i];
		container[i][13]+=(double)ss_n[2][i];
		container2[i][13]+=(double)ss_n[2][i]*ss_n[2][i];
	    }
    	}
    }else{
		if(opt_clust==1)
		{
	    	double** c=w_graph_compute_clust(WG, N_nodes);
	    	for(i=0;i<N_nodes;i++)
	    	{
	    		if(s[0][i] > 0)
	    		{		
					container[i][7]+=c[0][i];
					container2[i][7]+=c[0][i]*c[0][i];
					container[i][8]+=c[1][i];
					container2[i][8]+=c[1][i]*c[1][i];
	    		}
	    	}
	    	free_mat_double(c,2);
		}
		for(i=0;i<N_nodes;i++)
		{
			if(s[0][i] > 0)
	    	{
				//assert(s[0][i]>=k[0][i]);		
				node_nonzero[i][0]+=1;
				container[i][0]+=(double)k[0][i];
				container2[i][0]+=(double)k[0][i]*k[0][i];
				container[i][2]+=(double)s[0][i];
				container2[i][2]+=(double)s[0][i]*s[0][i];
				container[i][3]+=(double)y2[0][i];
				container2[i][3]+=(double)y2[0][i]*y2[0][i];
	    		container[i][4]+=(double)kk_n[0][i];
				container2[i][4]+=(double)kk_n[0][i]*kk_n[0][i];
				container[i][5]+=(double)kkw_n[0][i];
				container2[i][5]+=(double)kkw_n[0][i]*kkw_n[0][i];
				container[i][6]+=(double)ss_n[0][i];
				container2[i][6]+=(double)ss_n[0][i]*ss_n[0][i];
/*
		container[i][4]+=(double)xy[0][i]; // x
		container2[i][4]+=(double)xy[0][i]*(double)xy[0][i];
		container[i][5]+=(double)xy[1][i]; // y
		container2[i][5]+=(double)xy[1][i]*(double)xy[1][i];
		container[i][6]+=(double)xy[2][i]; //xy
		container2[i][6]+=(double)xy[2][i]*(double)xy[2][i];
*/
	    }
	}
	}
    T_container[0]+=T;
    T_container[1]+=T*T;
    free_mat_int(k,2);
    free_mat_int(s,2);
    free_mat_double(ss_n,4);
    free_mat_double(kk_n,4);
    free_mat_double(kkw_n,4);
    free_mat_double(y2,2);
    return;
}


void w_graph_node_stats_ensemble_print(int reps, int N_nodes, double* Tcont, double** cont, double ** cont2, int** node_nonzero, double av_k, double bin_exp, int len_acc, int opt_dir){
    char cadena[100];
    int i,j;
    //printf("i am printing\n");fflush(stdout);
    scale_vec_double(Tcont,1./reps,2);
    //average_matrix(cont, N_nodes, len_acc, reps);
    //average_matrix(cont2, N_nodes, len_acc, reps);
    ///// Only over existing but not for degrees or strengths!!!! ////
	//printf("alohaaa %d\n\n",len_acc/2);fflush(stdout);
    for(i=0;i<N_nodes;i++)
    {
	if (opt_dir==1)
	{
	    for(j=0;j<3;j++)
	    {
		cont[i][j]=cont[i][j]/reps;
		cont2[i][j]=cont2[i][j]/reps;
	    }	    
	    for(j=len_acc/2;j<len_acc/2+3;j++)
	    {
		cont[i][j]=cont[i][j]/reps;
		cont2[i][j]=cont2[i][j]/reps;
	    }	
	    
	    for(j=3;j<len_acc/2;j++)
	    {
			if(node_nonzero[i][0]>0)
			{
				cont[i][j]=cont[i][j]/node_nonzero[i][0];
				cont2[i][j]=cont2[i][j]/node_nonzero[i][0];
			}else{
				cont[i][j]=0;
				cont2[i][j]=0;
			}
	    }	    
	    for(j=len_acc/2+3;j<len_acc;j++)
	    {
			if(node_nonzero[i][1]>0)
			{
				cont[i][j]=cont[i][j]/node_nonzero[i][1];
				cont2[i][j]=cont2[i][j]/node_nonzero[i][1];
			}else{
				cont[i][j]=0;
				cont2[i][j]=0;
			}
	    }	
	}else{
	    for(j=0;j<3;j++)
	    {
		cont[i][j]=cont[i][j]/reps;
		cont2[i][j]=cont2[i][j]/reps;
	    }	    
	    for(j=3;j<len_acc;j++)
	    {
			if(node_nonzero[i][0]>0)
			{
				cont[i][j]=cont[i][j]/node_nonzero[i][0];
				cont2[i][j]=cont2[i][j]/node_nonzero[i][0];
			}else{
				cont[i][j]=0;
				cont2[i][j]=0;
			}
	    }
	}
    }
    if(opt_dir==1)
    {
		sprintf(cadena,"N%davs%.5f_ens_r%dnode_list.list",N_nodes,av_k,reps);
	}else{
		sprintf(cadena,"N%davs%.5f_undir_ens_r%d_node_list.list",N_nodes,av_k,reps);		
	}
    FILE* fil=open_file("w", cadena);
    fprintf(fil,"# <T>=%f+-%f # \n",Tcont[0] ,sqrt(Tcont[1]-Tcont[0]*Tcont[0]));
    if(opt_dir==1)
    {
    	fprintf(fil,"# Node_num\tk\tk_anal\tsset\tY2\tk_nn\tk^w_nn\ts^w_nn (out) then (in) # \n");
    }else{
	fprintf(fil,"# Node_num\tk\tk_anal\tsset\tY2\tk_nn\tk^w_nn\ts^w_nn (optionally \tc\tc^w) # \n");
    }
    //int i,j;
    double** node_atts=cast_mat_double(len_acc,N_nodes);
    for(i=0;i<N_nodes;i++)
    {
        fprintf(fil,"%d",i);
		//assert(cont[i][1]<=cont[i][2]);
		//assert(cont[i][7]<=cont[i][9]);
    	for(j=0;j<len_acc;j++)        
    	{
	    fprintf(fil," %f %f",cont[i][j],sqrt(cont2[i][j]-cont[i][j]*cont[i][j]));
	    node_atts[j][i] = cont[i][j];
	}
	    fprintf(fil,"\n");
    }
    fclose(fil);
    /*
    gsl_histogram * h1;
    // str in
    h1=histogram_double_log(WG->node_atts[2],0,max_value_double(WG->node_atts[2],N_nodes),bin_exp,N_nodes);
    sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%dnode_sIN.hist",N_nodes,av_k,expo-1,reps);
    print_acc(cadena, h1, h1);
    gsl_histogram_free(h1);
    // degrees in
    h1=histogram_double_log(WG->node_atts[0],0,max_value_double(WG->node_atts[0],N_nodes),bin_exp,N_nodes);
    sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%dnode_kIN.hist",N_nodes,av_k,expo-1,reps);
    print_acc(cadena, h1, h1);
    gsl_histogram_free(h1);

    if(opt_dir==1)
    {
    // degrees out
    h1=histogram_double_log(WG->node_atts[7],0,max_value_double(WG->node_atts[7],N_nodes),bin_exp,N_nodes);	
    sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%dnode_kOUT.hist",N_nodes,av_k,expo-1,reps);
    print_acc(cadena, h1, h1);
    gsl_histogram_free(h1);
    // str out
    h1=histogram_double_log(WG->node_atts[9],0,max_value_double(WG->node_atts[9],N_nodes),bin_exp,N_nodes);
    sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%dnode_sOUT.hist",N_nodes,av_k,expo-1,reps);
    print_acc(cadena, h1, h1);
    gsl_histogram_free(h1);
    }
    */
    return;
}
// P(k),P(s),P(w)

gsl_histogram ** w_graph_all_stats_ensemble_allocate(int dir, int s_min, int s_max, int k_min, int k_max, int w_max){
	int len_acc=2;
	int bins;
	bins=w_max;
	/*
	if(dir==1) 
	{
		len_acc=10;
	}else{
		len_acc=6;
	}
	*/
	gsl_histogram ** acc = (gsl_histogram**)malloc(sizeof(gsl_histogram*)*len_acc);
	acc[0] = set_acc_double(0,w_max,bins);
	acc[1] = set_acc_double(0,w_max,bins);
	/*
	if(dir==1)
	{
	
		acc[0] = set_acc_int(s_min,s_max);
		acc[1] = set_acc_int(s_min,s_max);
		acc[2] = set_acc_int(s_min,s_max);
		acc[3] = set_acc_int(s_min,s_max);

		acc[4] = set_acc_int(k_min,k_max);
		acc[5] = set_acc_int(k_min,k_max);
		acc[6] = set_acc_int(k_min,k_max);
		acc[7] = set_acc_int(k_min,k_max);

		acc[8] = set_acc_int(0,w_max);
		acc[9] = set_acc_int(0,w_max);

	}else{
		acc[0] = set_acc_int(s_min,s_max);
		acc[1] = set_acc_int(s_min,s_max);
		acc[2] = set_acc_int(k_min,k_max);
		acc[3] = set_acc_int(k_min,k_max);
		acc[4] = set_acc_int(0,w_max);
		acc[5] = set_acc_int(0,w_max);
	}
	*/
	return acc;
}

void w_graph_all_stats_ensemble_update(gsl_histogram** acc, W_GRAPH* WG, int N_nodes, int dir){
    	int E;
	//int **k=w_graph_compute_k(WG->node, N_nodes);
	//int **s=w_graph_compute_s(WG->node, N_nodes);
	int *w=w_graph_compute_w(WG, N_nodes, &E, -1);
	//long int T=sum_vec_int(s[0],N_nodes);
	/*
	if(dir==1)
	{
		update_int_acc(s[0], N_nodes, acc[0], acc[1], -1);
		update_int_acc(s[1], N_nodes, acc[2], acc[3], -1);
		update_int_acc(k[0], N_nodes, acc[4], acc[5], -1);
		update_int_acc(k[1], N_nodes, acc[6], acc[7], -1);
		update_int_acc(w, N_nodes, acc[8], acc[9], -1);
	}else{
		update_int_acc(s[0], N_nodes, acc[0], acc[1], -1);
		update_int_acc(k[0], N_nodes, acc[2], acc[3], -1);
		update_int_acc(w, N_nodes, acc[4], acc[5], -1);
	}
	*/
	update_int_acc(w, E, acc[0], acc[1], 1);
	free(w);
	return;
}
void w_graph_all_stats_ensemble_print(gsl_histogram** acc, int len, int reps, int N_nodes, double av_k, int opt_dir){
	char	cadena[100];
	//Normalize and std
	acc_normalize(acc,len, 1./(double)reps);
	acc_compute_std(acc,len, 1.);
	// Print all
	/*
	if(dir==1)
	{
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_sIN.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[0], acc[1]);
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_sOUT.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[2], acc[3]);
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_kIN.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[4], acc[5]);
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_kOUT.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[6], acc[7]);
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_w.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[8], acc[9]);
	}else{
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_sOUT.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[0], acc[1]);
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_kOUT.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[2], acc[3]);
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_w.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[4], acc[5]);
	}
	*/
    if(opt_dir==1)
    {
		sprintf(cadena,"N%davs%8.5f_ens_r%d_w.hist",N_nodes,av_k,reps);
	}else{
		sprintf(cadena,"N%davs%8.5f_undir_ens_r%d_w.hist",N_nodes,av_k,reps);
	}
	print_acc(cadena, acc[0], acc[1]);
	//Free all
	acc_free_all(acc, len);
	return;
}

