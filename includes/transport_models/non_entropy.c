/************************************************************
 *
 *                    Transport Models
 *
 *	Null models for weighted transport nets
 *
 *************************************************************/

#include "non_entropy.h"

/*********************************************/
/************** Lenormand et altr. (PLOS ONE) ********************/
/*
*    Lenormand M, Huet S, Gargiulo F, Deffuant G (2012) A Universal Model of
*    Commuting Networks. PLoS ONE 7(10): e45985. doi:10.1371/journal.pone.0045985
*/
/*********************************************/

double compute_gamma_frenchies(double surface){
	double alpha=3.15e-4;
	double nu=0.1777;
	double res;
	res=alpha*pow(surface,-nu);
	return res;
}
/**********************************/
double compute_p_french(int* s_in, double **d, double gamma, int dest, int origin, int* ind_in , int * ind_out, int len, int self_opt, int verbose){
	double norm=1e-15;
	double aux;
	int i,in_d,out_d,real_dest,real_origin;
	//compute numerator
	real_dest=ind_in[dest];
	real_origin=ind_out[origin];
	out_d=maxeq_int(real_dest,real_origin);
	in_d=mineq_int(real_dest,real_origin);
	double p=exp(-gamma*d[out_d][in_d])*(double)s_in[dest];
	for(i=0;i<len;i++)
	{
		real_dest=ind_in[i];
		if((real_dest!=real_origin) || (self_opt>0)) // avoid self-loops (or not)
		{
			out_d=maxeq_int(real_dest,real_origin);
			in_d=mineq_int(real_dest,real_origin);
			aux=(double)s_in[i]*exp(-d[out_d][in_d]*gamma);
			norm+=aux;
		}
	}
	if(verbose>0)
    {
        if(len<3) printf("Few to go, p/norm: %lf\n",p/norm);
        if(p>norm)
        {
            printf("Warning, p: %lf norm %lf Ins_av: %d\n",p,norm,len); fflush(stdout);
        }
    }
	p/=norm;
	return p;
}
/**********************************/
void check_norm_french(int *s_in, double **d, double gamma, int origin, int* ind_in, int *ind_out, int len, int self_opt, int verbose){
	double norm=0.;
	int real_dest,i;
	int real_origin=ind_out[origin];
	//int out_d,in_d;
	for(i=0;i<len;i++)
	{
		real_dest=ind_in[i];
		if((real_dest!=real_origin) || (self_opt>0)) // accept self-loops or not
        {
            norm+=compute_p_french(s_in, d, gamma, i , origin, ind_in, ind_out, len);
        }
	}
	if(fabs(norm-1.)>=1e-10)
	{
        if(verbose>0)
        {
            printf("Normalization is not right (%lf) ! Aborting...\n",norm);
            printf("Info: Node origin %i, length of sum %i \n",real_origin,len);
        }
		/*for(i=0;i<len;i++)
		{
			out_d=maxeq_int(ind_in[i],real_origin);
			in_d=mineq_int(ind_in[i],real_origin);
			printf("s_in : %d d_ij : %lf p: %lf gamma: %lf\n",s_in[i],d[out_d][in_d],compute_p(s_in, d, gamma, i , origin, ind_in, ind_out, len),gamma);
		}*/
		abort();
	}
	return;
}
/**********************************/
//// i am here !!! ///////
w_graph* w_graph_seq_gravity_directed_graph(int N_nodes,int **s, double **d, double gamma, gsl_rng* randgsl, int Max_fails, int self_opt, int verbose){
    // alloc w_graph
    w_graph* node = malloc(sizeof(w_graph)*N_nodes);
    aux=0;
    for(i=0;i<N_nodes;i++)
    {
        node[i].idnum=i;
        node[i].kin=0;
        node[i].kout=0;
        node[i].kout=0;
        node[i].sin=0;
        node[i].sout=0;
        node[i].mem_in=1;
        node[i].mem_out=1;
        node[i].out=calloc(1,sizeof(int));
        node[i].in=calloc(1,sizeof(int));
        node[i].w_in=calloc(1,sizeof(int));
        node[i].w_out=calloc(1,sizeof(int));
    }
    // allocs
	int* s_in=cast_vec_int(N_nodes);
	int* s_inf=cast_vec_int(N_nodes);
	int* s_out=cast_vec_int(N_nodes);
	int* s_outf=cast_vec_int(N_nodes);
	int* in_inds=cast_vec_int(N_nodes);
	int* out_inds=cast_vec_int(N_nodes);
    
	int origin,dest,aux1,aux2,s_fake,s_fake2;
	int out_av,in_av,aux_i,aux_j;
    int l,l2;
    
	double ind1,ind2,ind3;
	ind1=ind2=ind3=0;

	double ind11,ind22,ind33;
	ind11=ind22=ind33=0;

	double p,r1;

	int tot_trips,tot_trips2,tot_trips3,flag;
 	int* in_finds=cast_vec_int(N_nodes);
	int* out_finds=cast_vec_int(N_nodes);

    in_av=aux1=0;
	out_av=aux2=0;
	tot_trips2=tot_trips3=0;
	
    // check all good
    for(i=0;i<N_nodes;i++)
	{
		if(s[1][i]!=0)
		{
			s_inf[aux1]=s[1][i];
			in_finds[i]=aux1;
			aux1++;
			//tot_trips2+=s[i][1];
		}else{
			s_inf[N_nodes-in_av-1]=0;
			in_finds[i]=N_nodes-1-in_av;
			in_av++;
		}
		if(s[0][i]!=0)
		{
			s_outf[aux2]=s[0][i];
			out_finds[i]=aux2;
			aux2++;
			//tot_trips3+=s[i][0];
		}else{
			s_outf[N_nodes-out_av-1]=0;
			out_finds[i]=N_nodes-1-out_av;
			out_av++;
		}
	}
	tot_trips2=sum_vec_int(s_outf,N_nodes);
	tot_trips3=sum_vec_int(s_inf,N_nodes);
	for(i=0;i<N_nodes;i++)
	{
		l=count_appearances_int(in_finds,i,N_nodes);
		l2=count_appearances_int(out_finds,i,N_nodes);
		if((l!=1) || (l2!=1)) // check 1 to 1 mapping
		{
			if(verbose>0)printf("Incorrect mapping for i: %i %i %i\n",i,l,l2);fflush(stdout);
			abort();
		}
	}
	if(verbose>0)
    {
        if(tot_trips2-tot_trips3 != 0)printf("Warning, total outs (%i) and total ins (%i) different!\n",tot_trips2,tot_trips3);
        printf("Available in %i, Available out %i\n",N_nodes-in_av,N_nodes-out_av);fflush(stdout);
    }
	//for(i=0;i<N_nodes;i++){
	//	if(s_outf[i]==0 || s_inf[i]==0){printf("%i %i %i %i\n",(int)out_finds[i],(int)s_outf[i],(int)in_finds[i],(int)s_inf[i]);fflush(stdout);}
	//}


	//go!
	int fails=0;
	int good_reps=0;
	int in_av_def,out_av_def;
	in_av_def=N_nodes-in_av;
	out_av_def=N_nodes-out_av;
	for(r=0;r<Reps;r++){
		flag=-1;
		tot_trips2=0;
		// Initial conditions
		for(i=0;i<N_nodes;i++)
		{
			s_in[i]=s_inf[i]; // copy ins
			s_out[i]=s_outf[i]; // copy outs
			in_inds[i]=in_finds[i]; // copy indices
			out_inds[i]=out_finds[i]; // copy indices
			for(j=0;j<N_nodes;j++)
			{ // set row to 0
				w_false[i][j]=0;
			}
		}
		//out_av=in_av=500;
		out_av=out_av_def;
		in_av=in_av_def;
		//in_av=N_nodes-count_appearances_int(s_in,0,N_nodes);
		//out_av=N_nodes-count_appearances_int(s_out,0,N_nodes);
		do{// While available outs
			aux1=aux2=-1;
			do
			{
				origin  = (int)gsl_rng_uniform_int(randgsl,(unsigned int)out_av); //out from out_av
				dest    = (int)gsl_rng_uniform_int(randgsl,(unsigned int)in_av); // in from in_av
				aux1=out_inds[origin];
				aux2=in_inds[dest];
			}while(aux1==aux2);
				//p  = compute_p_french(s_in, d, gamma, dest, origin, in_inds, out_inds, in_av);
				p = 1./(double)in_av;
			/*if(in_av>1)
			{
				p  = compute_p_french(s_in, d, gamma, dest, origin, in_inds, out_inds, in_av);
			}
			else{
				p=1.;
			}*/
			//check_norm_french(s_in,d,gamma,origin, in_inds,out_inds,in_av);
			r1      = gsl_rng_uniform(randgsl);
			if(r1<=p)
			{
				fails=0;
				i=out_inds[origin];
				j=in_inds[dest];
				w_false[i][j]++;
				tot_trips2++;
				s_in[dest]--;
				s_out[origin]--;
				//printf("Accepted trip! From %i to %i, s_out %i s_in %i To go %i\n",i,j,s_out[origin],s_in[dest],tot_trips3-tot_trips2);fflush(stdout);
				if(s_in[dest]<1)
				{
					//printf("%i",s_in[dest]);
					s_fake=s_in[in_av-1];
					s_fake2=s_in[dest];
					s_in[dest]=s_fake;
					s_in[in_av-1]=s_fake2;
					//printf(" %i \n",s_in[dest]);fflush(stdout);
					aux_i=in_inds[dest];
					aux_j=in_inds[in_av-1];
					in_inds[in_av-1]=aux_i;
					in_inds[dest]=aux_j;
					in_av--;
					//printf("%d Exhausted. Available in %d, Available out %d p %lf 2go: %d\n",j,in_av,out_av,p,tot_trips3-tot_trips2);fflush(stdout);
					//printf("%i %d\n",count_appearances_int(s_in,0,N_nodes),count_appearances_int(s_out,0,N_nodes));fflush(stdout);
					//printf(" ++ Sink exhausted! To go: %i \n",tot_trips3-tot_trips2);fflush(stdout);
				}
				if(s_out[origin]<1)
				{
					//printf("%i",s_out[origin]);
					s_fake=s_out[out_av-1];
					s_fake2=s_out[origin];
					s_out[origin]=s_fake;
					s_out[out_av-1]=s_fake2;
					//printf(" %i \n",s_out[origin]);fflush(stdout);
					aux_i=out_inds[origin];
					aux_j=out_inds[out_av-1];
					out_inds[out_av-1]=aux_i;
					out_inds[origin]=aux_j;
					out_av--;
					//printf("%d Exhausted. Available in %i, Available out %i p %lf 2go: %d\n",i,in_av,out_av,p,tot_trips3-tot_trips2);fflush(stdout);
					//printf("%i %i\n",count_appearances_int(s_in,0,N_nodes),count_appearances_int(s_out,0,N_nodes));fflush(stdout);
					//printf(" -- Source exhausted! To go: %i \n",tot_trips3-tot_trips2);fflush(stdout);
				}
			}
			else
			{
				fails++;
				//printf("Max_fails %i\n",fails);
				//fflush(stdout);
			}
			if((out_av<=0) || (in_av<=0))
			{
				flag=0;
			}
			if((out_av==1) && (in_av==1) && (out_inds[0]==in_inds[0]))
			{
				flag=1;// avoid last self loop
				printf("Last trips are self-loops, ignoring last %i trips!\n",tot_trips3-tot_trips2);
			}
			if(tot_trips3-tot_trips2<0)
			{
				flag=3;
				printf("Something is wrong: too many trips!\n");
				abort();
			}
			if(fails>Max_fails)
			{
				flag=2;// fail
				printf("Maximum rejections attained, ignoring!\n");
			}
		}while(flag<0);
		if(flag!=2)
		{
			good_reps++;
			// indices
			p=sorensen_index(w_real,w_false,N_nodes);
			ind1+=p;
			ind11+=p*p;
			p=nmae_index(w_real,w_false,N_nodes);
			ind2+=p;ind22+=p*p;
			p=nrmse_index(w_real,w_false,N_nodes);
			ind3+=p;ind33+=p*p;
			// add to accumulator
			tot_trips=0;
			for(i=0;i<N_nodes;i++)
			{
				for(j=0;j<N_nodes;j++)
				{
					if(w_false[i][j]>0)
					{
						w[i][j]=w[i][j]+(double)w_false[i][j];
						w2[i][j]=w2[i][j]+(double)(w_false[i][j]*w_false[i][j]);
						tot_trips+=w_false[i][j];
					}
				}
			}
			printf("Total trips for this round: %i, Suposed to be: %i\n",tot_trips,tot_trips3);
			rep_count++;
			if(rep_count%100 == 0) {printf("*** Done %i of %i reps\n **** ",rep_count,Reps);fflush(stdout);}
			if(rep_count == 1) 
			{
				printf("*** Printing sample file\n **** ");fflush(stdout);
				char* address="single_sample_french.sample";
				print_w_2_file_int(address,N_nodes,w_false);
			}
		}
	}
	ind1=ind1/(double)good_reps;
	ind11=sqrt(ind11/(double)good_reps-ind1*ind1);
	ind2=ind2/(double)good_reps;
	ind22=sqrt(ind22/(double)good_reps-ind2*ind2);
	ind3=ind3/(double)good_reps;
	ind33=sqrt(ind33/(double)good_reps-ind3*ind3);
	write_indices(ind1,ind11,ind2,ind22,ind3,ind33,output_indices);
	// get average
	for(i=0;i<N_nodes;i++)
	{
		free(w_false[i]);
		for(j=0;j<N_nodes;j++)
		{
			w[i][j]=w[i][j]/(double)Reps;
			w2[i][j]=sqrt(w2[i][j]/(double)Reps-w[i][j]*w[i][j]);
			//printf("%lf %lf %i %i\n",w[i][j],w2[i][j],i,j);fflush(stdout);
		}
	}
	print_w_2_file_double(output_name,N_nodes,w,w2);
	free(s_in);free(s_out);free(in_inds);free(out_inds);
	free(s_inf);free(s_outf);free(in_finds);free(out_finds);
	free(w_false);
	for(i=0;i<N_nodes;i++)
	{
		free(w[i]);
		free(w2[i]);
	}
	free(w);
	free(w2);
	return;
} 
/*********************************************/
/************** Radiation ********************/
/*********************************************/
int destination_rad_model(int** s, double** d, int origin, int N_nodes, int T, gsl_rng * randgsl)
{
	int i,j,dest;
	double max_val,r,node_val,dist2,dist=100000000.,eps=1e-15;
	max_val=0.;
	int finite=rint((double)s[origin][0]/(1.-((double)s[origin][1]/(double)T)));
	for(j=0;j<finite;j++)
	{//generate randoms for original node
		r=gsl_rng_uniform(randgsl);
		if(r>max_val)max_val=r;
	}
	node_val=max_val;
	dest=origin;
	for(i=0;i<N_nodes;i++)
	{//for each node (others)
		if(i!=origin)
		{
			max_val=0;
			for(j=0;j<s[i][1];j++)
			{//generate s_in randoms
				r=gsl_rng_uniform(randgsl);
				if(r>max_val)max_val=r;
			}
			if(max_val>=node_val)
			{//keep the ones that exceed threshold
				dist2=d[maxeq_int(origin,i)][mineq_int(origin,i)]; //store dist
				if((dist2+eps)<=dist) // if closer
				{
					if(fabs(dist-dist2)<eps)//if equal undraw
					{
						r=gsl_rng_uniform(randgsl);
						if(r>0.5)
						{
							dist=dist2;
							dest=i;
						}
					}else{
						dist=dist2;
						dest=i;
					}
				}
			}
		}
	}
	return dest;
}
/*********************************************/
void radiation_model(int **w_real,int N_nodes,int **s, double **d, int Reps,int seed, char *output_name, char *output_indices)
{
	gsl_rng * randgsl = gsl_rng_alloc(gsl_rng_taus);	/// we initialize the random generator of the gsl
	gsl_rng_set(randgsl,seed);
	double** w=(double**)malloc(sizeof(double*)*N_nodes);//w_matrix w[i][j][0:1]
	double** w2=(double**)malloc(sizeof(double*)*N_nodes);//w_matrix w[i][j][0:1]
	int** w_false=(int**)malloc(sizeof(int*)*N_nodes);//w_fake matrix w[i][j]
	int i,j,r; // aux indices
	for(i=0;i<N_nodes;++i)
	{
		w[i]=(double*)malloc(sizeof(double*)*(N_nodes)); // column
		w2[i]=(double*)malloc(sizeof(double*)*(N_nodes)); // column
		w_false[i]=(int*)malloc(sizeof(int)*(N_nodes)); // full matrix
		for(j=0;j<N_nodes;++j)
		{
			w[i][j]=0.;
			w2[i][j]=0.;
			w_false[i][j]=0;
		
		}
	}
	// Go for reps
	double ind1,ind2,ind3,ind11,ind22,ind33,aux;
	ind1=ind2=ind3=ind11=ind22=ind33=0;
	int tot_trips2,tot_trips3;
	int given_trips,dest;
	tot_trips2=tot_trips3=0;
	for(i=0;i<N_nodes;i++)
	{
		tot_trips2+=s[i][0];
		tot_trips3+=s[i][1];
	}
	if(tot_trips2-tot_trips3 != 0)
	{printf("Warning, total outs (%i) and total ins (%i) different!\n",tot_trips2,tot_trips3);}
	for(r=0;r<Reps;r++)
	{
		tot_trips2=0;
		// Initial conditions
		for(i=0;i<N_nodes;i++)
		{
			for(j=0;j<N_nodes;j++)w_false[i][j]=0;
		}
		for(i=0;i<N_nodes;i++) // for each node (they are indep)
		{
			given_trips=s[i][0];
			//printf("Doing %i %i\n",i,given_trips);fflush(stdout);
			while(given_trips>0) // if it has indeed out-going trips!
			{
				dest = 	destination_rad_model(s, d, i, N_nodes,tot_trips3,  randgsl);
				given_trips--; //regardless of success
				if(dest!=i) //if success
				{
					w_false[i][dest]+=1;
					tot_trips2+=1;
				}
			}
			if(i%100==0)printf("---Performed %i nodes out of %i, %i trips to go---\n",i,N_nodes,tot_trips3-tot_trips2);
		}
		// indices
		aux=sorensen_index(w_real,w_false,N_nodes);
		ind1+=aux;ind11+=aux*aux;
		aux=pnzt_index(w_real,w_false,N_nodes);
		ind2+=aux;ind22+=aux*aux;
		aux=nrmse_index(w_real,w_false,N_nodes);
		ind3+=aux;ind33+=aux*aux;
		// add to accumulator
		tot_trips2=0;
		for(i=0;i<N_nodes;i++)
		{
			for(j=0;j<N_nodes;j++)
			{
				w[i][j]=w[i][j]+(double)w_false[i][j];
				w2[i][j]=w2[i][j]+(double)(w_false[i][j]*w_false[i][j]);
				tot_trips2+=w_false[i][j];
			}
		}
		printf("Total trips for this round: %i, Suposed to be: %i\n",tot_trips2,tot_trips3);
		if(r%100 == 0)
		{
			printf("*** Done %i of %i reps\n **** ",r,Reps);fflush(stdout);
		}
		if(r == 1) 
		{
			printf("*** Printing sample file\n **** ");fflush(stdout);
			char* address="single_sample_rad_stoc.sample";
			print_w_2_file_int(address,N_nodes,w_false);
		}
	}
	ind1=ind1/(double)Reps;
	ind11=sqrt(ind11/(double)Reps-ind1*ind1);
	ind2=ind2/(double)Reps;
	ind22=sqrt(ind22/(double)Reps-ind2*ind2);
	ind3=ind3/(double)Reps;
	ind33=sqrt(ind33/(double)Reps-ind3*ind3);
	write_indices(ind1,ind11,ind2,ind22,ind3,ind33,output_indices);
	// get average
	for(i=0;i<N_nodes;i++)
	{
		free(w_false[i]);
		for(j=0;j<N_nodes;j++)
		{
			w[i][j]=w[i][j]/(double)Reps;
			w2[i][j]=sqrt(w2[i][j]/(double)Reps-w[i][j]*w[i][j]);
		}
	}
	free(w_false);
	print_w_2_file_double(output_name,N_nodes,w,w2);
	for(i=0;i<N_nodes;i++)
	{
		free(w[i]);
		free(w2[i]);
		free(s[i]);
	}
	free(s);
	free(w);
	free(w2);
	return;
}
/*********************************************/
void radiation_model_multinomial(int **w_real,int N_nodes,int **s, double **d, int Reps,int seed, char *output_name, char *output_indices)
{
	gsl_rng * randgsl = gsl_rng_alloc(gsl_rng_taus);	/// we initialize the random generator of the gsl
	gsl_rng_set(randgsl,seed);
	double** w=(double**)malloc(sizeof(double*)*N_nodes);//w_matrix w[i][j][0:1]
	double** w2=(double**)malloc(sizeof(double*)*N_nodes);//w_matrix w[i][j][0:1]
	int** w_false=(int**)malloc(sizeof(int*)*N_nodes);//w_fake matrix w[i][j]
	int i,j,r,dummy_count; // aux indices
	for(i=0;i<N_nodes;++i)
	{
		w[i]=(double*)malloc(sizeof(double*)*(N_nodes)); // column
		w2[i]=(double*)malloc(sizeof(double*)*(N_nodes)); // column
		w_false[i]=(int*)malloc(sizeof(int)*(N_nodes)); // full matrix
		for(j=0;j<N_nodes;++j)
		{
			w[i][j]=0.;
			w2[i][j]=0.;
			w_false[i][j]=0;
		}
	}
	// Go for reps
	double ind1,ind2,ind3,ind11,ind22,ind33,aux;
	ind1=ind2=ind3=ind11=ind22=ind33=0;
	int tot_trips2,tot_trips3,sij;
	tot_trips2=tot_trips3=0;
	double** ps=cast_mat_double(N_nodes,N_nodes);
	unsigned int** ps_false=(unsigned int**)malloc(sizeof(unsigned int*)*N_nodes);
	int aux_int=0;
	for(i=0;i<N_nodes;i++)
	{
		ps_false[i]=malloc(sizeof(unsigned int)*N_nodes);
		tot_trips2+=s[i][0];
		tot_trips3+=s[i][1];
	}

	if(tot_trips2-tot_trips3 != 0)
	{printf("Warning, total outs (%i) and total ins (%i) different!\n",tot_trips2,tot_trips3);}
	for(i=0;i<N_nodes;i++) // each node
	{
		for(j=0;j<N_nodes;j++) // each available node
		{
			if((s[i][1]>0) && (s[j][1]>0) && (s[i][0]>0) &&(i!=j))
			{
				sij=compute_sij_rad(s,d,i,j,N_nodes);
				ps[i][j]=(double)(s[i][1]*s[j][1])/((double)(s[i][1]+sij)*(double)(s[i][1]+sij+s[j][1])); // probabilities
				ps[i][j]=ps[i][j]*1./(1.-((double)s[i][1]/(double)tot_trips2));
				//ps[i][j]=ps[i][j]*(double)s[i][0];
			}
			else
			{
				ps[i][j]=0;
			}
			ps_false[i][j]=0.;
		}
		ps_false[i][i]=0;
		ps[i][i]=0.;
	}
	for(r=0;r<Reps;r++)
	{
		tot_trips2=0;
		// Initial conditions
		//go!
		for(i=0;i<N_nodes;i++)
		{
			dummy_count=0;
			w_false[i][i]=0;
			gsl_ran_multinomial(randgsl, N_nodes, s[i][0], ps[i], ps_false[i]); // run multinomial
			for(j=0;j<N_nodes;j++)
			{
				if(i!=j)
				{
				//if(ps_false[aux_int]>0) printf("%d %d %d %f\n",i,j,ps_false[aux_int],ps[aux_int]);
				//w_false[i][j]=gsl_ran_binomial ( randgsl, ps[i][j], s[i][0]);
				w_false[i][j]=ps_false[i][j];
				aux_int++;
				dummy_count+=w_false[i][j];
				}
			}
			//printf("Allocated trips for i:%i Talloc_i: %i  T_i:%i\n",i,dummy_count-w_false[i][i],s[i][0]);fflush(stdout);
		}
		// indices
		aux=sorensen_index(w_real,w_false,N_nodes);
		ind1+=aux;ind11+=aux*aux;
		aux=pnzt_index(w_real,w_false,N_nodes);
		ind2+=aux;ind22+=aux*aux;
		aux=nrmse_index(w_real,w_false,N_nodes);
		ind3+=aux;ind33+=aux*aux;
		// add to accumulator
		tot_trips2=0;
		for(i=0;i<N_nodes;i++)
		{
			for(j=0;j<N_nodes;j++)
			{
				w[i][j]=w[i][j]+(double)w_false[i][j];
				w2[i][j]=w2[i][j]+(double)(w_false[i][j]*w_false[i][j]);
				tot_trips2+=w_false[i][j];
			}
		}
		printf("Total trips for this round: %i, Suposed to be: %i\n",tot_trips2,tot_trips3);
		if(r%100 == 0)
		{
			printf("*** Done %i of %i reps\n **** ",r,Reps);fflush(stdout);
		}
		if(r == 1) 
		{
			printf("*** Printing sample file\n **** ");fflush(stdout);
			char* address="single_sample_rad_multi.sample";
			print_w_2_file_int(address,N_nodes,w_false);
		}
	}
	ind1=ind1/(double)Reps;
	ind11=sqrt(ind11/(double)Reps-ind1*ind1);
	ind2=ind2/(double)Reps;
	ind22=sqrt(ind22/(double)Reps-ind2*ind2);
	ind3=ind3/(double)Reps;
	ind33=sqrt(ind33/(double)Reps-ind3*ind3);
	write_indices(ind1,ind11,ind2,ind22,ind3,ind33,output_indices);
	// get average
	for(i=0;i<N_nodes;i++)
	{
		free(w_false[i]);
		for(j=0;j<N_nodes;j++)
		{
			w[i][j]=w[i][j]/(double)Reps;
			w2[i][j]=sqrt(w2[i][j]/(double)Reps-w[i][j]*w[i][j]);
		}
	}
	free(w_false);
	print_w_2_file_double(output_name,N_nodes,w,w2);
	for(i=0;i<N_nodes;i++)
	{
		free(w[i]);
		free(w2[i]);
		free(s[i]);
	}
	free(s);
	free(w);
	free(w2);
	return;
}
/*********************************************/
int compute_sij_rad(int  **s, double ** dist, int origin, int dest, int N_nodes){
	int i;
	int sij=0;
	double dd=dist[maxeq_int(origin,dest)][mineq_int(origin,dest)];
	for(i=0;i<N_nodes;i++)
	{
		if((i!=origin)&&(i!=dest)&&(dist[maxeq_int(i,origin)][mineq_int(i,dest)]<=dd))
		{
			sij+=s[i][1];
		}
	}
	return sij;
}



/*********************************************/
/*********** Bose Models  ********************/
/*********************************************/
void bose_model_dist(int **w_real,int N_nodes,int **s, double **x, double **d, double gamma, int Reps,int seed, char *output_name, char *output_indices)
{
	gsl_rng * randgsl = gsl_rng_alloc(gsl_rng_taus);	/// we initialize the random generator of the gsl
	gsl_rng_set(randgsl,seed);
	double** w=(double**)malloc(sizeof(double*)*N_nodes);//w_matrix w[i][j][0:1]
	double** w2=(double**)malloc(sizeof(double*)*N_nodes);//w_matrix w[i][j][0:1]
	int** w_false=(int**)malloc(sizeof(int*)*N_nodes);//w_fake matrix w[i][j]
	int i,j,r; // aux indices
	for(i=0;i<N_nodes;++i)
	{
		w[i]=(double*)malloc(sizeof(double*)*(N_nodes)); // column
		w2[i]=(double*)malloc(sizeof(double*)*(N_nodes)); // column
		w_false[i]=(int*)malloc(sizeof(int)*(N_nodes)); // full matrix
		for(j=0;j<N_nodes;++j)
		{
			w[i][j]=0.;
			w2[i][j]=0.;
			w_false[i][j]=0;
			
		}
	}
	// Go for reps
	double ind1,ind2,ind3,ind11,ind22,ind33,aux;
	double p;
	ind1=ind2=ind3=ind11=ind22=ind33=0;
	int tot_trips2,tot_trips3;
	int flag,s_node;
	tot_trips2=tot_trips3=0;
	for(i=0;i<N_nodes;i++)
	{
		tot_trips2+=s[i][0];
		tot_trips3+=s[i][1];
	}
	if(tot_trips2-tot_trips3 != 0)
	{printf("Warning, total outs (%i) and total ins (%i) different!\n",tot_trips2,tot_trips3);}
	for(r=0;r<Reps;r++)
	{
		tot_trips2=0;
		// Initial conditions
		for(i=0;i<N_nodes;i++)
		{
			for(j=0;j<N_nodes;j++)w_false[i][j]=0;
		}
		//go!
		for(i=0;i<N_nodes;i++) // for each node (they are indep)
		{	
			s_node=0;
			for(j=0;j<N_nodes;j++)
			{
				if(j!=i) // no self loops
				{
					flag=0;
					while(flag>0) // while not a fail
					{
						p=x[i][0]*x[j][1]*exp(-gamma*d[maxeq_int(i,j)][mineq_int(i,j)]);
						r=gsl_rng_uniform(randgsl);
						if(r<=p)
						{
							w_false[i][j]+=1;
							tot_trips2+=1;
							s_node++;
						}
						else
						{
							flag=-1;	
						}
					}
				}
			}
			printf("Node: %i Original Strength out:%i Strength out: %i Difference: %i",i,s[i][0],s_node,s[i][0]-s_node);
			if(i%100==0)printf("---Performed %i nodes out of %i, %i trips to go---\n",i,N_nodes,tot_trips3);
		}
		// indices
		aux=sorensen_index(w_real,w_false,N_nodes);
		ind1+=aux;ind11+=aux*aux;
		aux=pnzt_index(w_real,w_false,N_nodes);
		ind2+=aux;ind22+=aux*aux;
		aux=nrmse_index(w_real,w_false,N_nodes);
		ind3+=aux;ind33+=aux*aux;
		// add to accumulator
		tot_trips2=0;
		for(i=0;i<N_nodes;i++)
		{
			for(j=0;j<N_nodes;j++)
			{
				w[i][j]=w[i][j]+(double)w_false[i][j];
				w2[i][j]=w2[i][j]+(double)(w_false[i][j]*w_false[i][j]);
				tot_trips2+=w_false[i][j];
			}
		}
		printf("Total trips for this round: %i, Suposed to be: %i\n",tot_trips2,tot_trips3);
		if(r%100 == 0)
		{
			printf("*** Done %i of %i reps\n **** ",r,Reps);fflush(stdout);
		}
		if(r == 0) 
		{
			printf("*** Printing sample file\n **** ");fflush(stdout);
			char* address="single_sample_bose.sample";
			print_w_2_file_int(address,N_nodes,w_false);		}
	}
	ind1=ind1/(double)Reps;
	ind11=sqrt(ind11/(double)Reps-ind1*ind1);
	ind2=ind2/(double)Reps;
	ind22=sqrt(ind22/(double)Reps-ind2*ind2);
	ind3=ind3/(double)Reps;
	ind33=sqrt(ind33/(double)Reps-ind3*ind3);
	write_indices(ind1,ind11,ind2,ind22,ind3,ind33,output_indices);
	// get average
	for(i=0;i<N_nodes;i++)
	{
		free(w_false[i]);
		for(j=0;j<N_nodes;j++)
		{
			w[i][j]=w[i][j]/(double)Reps;
			w2[i][j]=sqrt(w2[i][j]/(double)Reps-w[i][j]*w[i][j]);
		}
	}
	free(w_false);
	print_w_2_file_double(output_name,N_nodes,w,w2);
	for(i=0;i<N_nodes;i++)
	{
		free(w[i]);
		free(w2[i]);
		free(s[i]);
	}
	free(s);
	free(w);
	free(w2);
	return;
}
/*********************************************/
/*********** Wilson Models  ********************/
/*********************************************/
/*********************************************/
void wilson_model_distance(int **w_real,int N_nodes,int **s, double **x, double **d, double gamma, int Reps,int seed, char *output_name, char *output_indices)
{
	gsl_rng * randgsl = gsl_rng_alloc(gsl_rng_taus);	/// we initialize the random generator of the gsl
	gsl_rng_set(randgsl,seed);
	double** w=(double**)malloc(sizeof(double*)*N_nodes);//w_matrix w[i][j][0:1]
	double** w2=(double**)malloc(sizeof(double*)*N_nodes);//w_matrix w[i][j][0:1]
	int** w_false=(int**)malloc(sizeof(int*)*N_nodes);//w_fake matrix w[i][j]
	int i,j,r; // aux indices
	for(i=0;i<N_nodes;++i)
	{
		w[i]=(double*)malloc(sizeof(double*)*(N_nodes)); // column
		w2[i]=(double*)malloc(sizeof(double*)*(N_nodes)); // column
		w_false[i]=(int*)malloc(sizeof(int)*(N_nodes)); // full matrix
		for(j=0;j<N_nodes;++j)
		{
			w[i][j]=0.;
			w2[i][j]=0.;
			w_false[i][j]=0;			
		}
	}
	// Go for reps
	double ind1,ind2,ind3,ind11,ind22,ind33,aux;
	ind1=ind2=ind3=ind11=ind22=ind33=0;
	int tot_trips2,tot_trips3;
	tot_trips2=tot_trips3=0;
	double* ps=(double*)malloc(sizeof(double)*N_nodes*(N_nodes-1));
	unsigned int* ps_false=(unsigned int*)malloc(sizeof(unsigned int)*N_nodes*(N_nodes-1));
	int aux_int=0;
	for(i=0;i<N_nodes;i++)
	{
		tot_trips2+=s[i][0];
		tot_trips3+=s[i][1];
	}

	if(tot_trips2-tot_trips3 != 0)
	{printf("Warning, total outs (%i) and total ins (%i) different!\n",tot_trips2,tot_trips3);}
	for(i=0;i<N_nodes;i++)
	{
		for(j=0;j<N_nodes;j++)
		{
			if(j!=i)
			{
				ps[aux_int]=x[i][0]*x[j][1]*exp(-gamma*d[maxeq_int(i,j)][mineq_int(i,j)]);
				ps_false[aux_int]=0.;
				aux_int++;
			}
		}
	}

	for(r=0;r<Reps;r++)
	{
		tot_trips2=0;
		// Initial conditions
		//go!
		gsl_ran_multinomial(randgsl, N_nodes*(N_nodes-1),  tot_trips3, ps, ps_false);
		aux_int=0;
		for(i=0;i<N_nodes;i++)
		{
			w_false[i][i]=0;
			for(j=0;j<N_nodes;j++)
			{
					if(j!=i)
				{
					//if(ps_false[aux_int]>0) printf("%d %d %d %f\n",i,j,ps_false[aux_int],ps[aux_int]);
					w_false[i][j]=ps_false[aux_int];
					aux_int++;
				}
			}
		}
		// indices
		aux=sorensen_index(w_real,w_false,N_nodes);
		ind1+=aux;ind11+=aux*aux;
		aux=pnzt_index(w_real,w_false,N_nodes);
		ind2+=aux;ind22+=aux*aux;
		aux=nrmse_index(w_real,w_false,N_nodes);
		ind3+=aux;ind33+=aux*aux;
		// add to accumulator
		tot_trips2=0;
		for(i=0;i<N_nodes;i++)
		{
			for(j=0;j<N_nodes;j++)
			{
				w[i][j]=w[i][j]+(double)w_false[i][j];
				w2[i][j]=w2[i][j]+(double)(w_false[i][j]*w_false[i][j]);
				tot_trips2+=w_false[i][j];
			}
		}
		printf("Total trips for this round: %i, Suposed to be: %i\n",tot_trips2,tot_trips3);
		if(r%100 == 0)
		{
			printf("*** Done %i of %i reps\n **** ",r,Reps);fflush(stdout);
		}
		if(r == 1) 
		{
			printf("*** Printing sample file\n **** ");fflush(stdout);
			char* address="single_sample_wilson.sample";
			print_w_2_file_int(address,N_nodes,w_false);
		}
	}
	ind1=ind1/(double)Reps;
	ind11=sqrt(ind11/(double)Reps-ind1*ind1);
	ind2=ind2/(double)Reps;
	ind22=sqrt(ind22/(double)Reps-ind2*ind2);
	ind3=ind3/(double)Reps;
	ind33=sqrt(ind33/(double)Reps-ind3*ind3);
	write_indices(ind1,ind11,ind2,ind22,ind3,ind33,output_indices);
	// get average
	for(i=0;i<N_nodes;i++)
	{
		free(w_false[i]);
		for(j=0;j<N_nodes;j++)
		{
			w[i][j]=w[i][j]/(double)Reps;
			w2[i][j]=sqrt(w2[i][j]/(double)Reps-w[i][j]*w[i][j]);
		}
	}
	free(w_false);
	print_w_2_file_double(output_name,N_nodes,w,w2);
	for(i=0;i<N_nodes;i++)
	{
		free(w[i]);
		free(w2[i]);
		free(s[i]);
	}
	free(s);
	free(w);
	free(w2);
	return;
}
/*********************************************/
// Check results for bose and for wilson!
