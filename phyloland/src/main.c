#include <stdio.h>
#include <stdlib.h>
#include "R.h"
#include <Rmath.h>


int IsFiniteNumber(double x) ;
long sample_proba( int x[], double proba[], long length);
double rexp_proba(double lambda_exp);
double dexp_proba(double x, double lambda_exp);
double distkm(double lat1, double lat2, double long1, double long2);

void migC(double *space, double *mat_Dists1, double *mat_Dists2, int *occupied, double *sigma, double *lambda , double *tau, double *mig, double *mig_event, int *space_dim, int *space_size, int *length_mig)
{
	int i=0, j=0, index=0;

	// Affichage des parametres //

	/*printf("\nspace_dim\n %d\n",*space_dim);
	printf("\nspace_size\n %d\n",*space_size);
	printf("\nspace\n");
	for ( i = 0; i < ((*space_size)*(*space_dim)) ; i++) {
		printf("%lf \t",space[i]);
	} */
	
	
	/*printf("\nmat_Dists1\n");
	 for ( i = 0; i < ((*space_size)) ; i++) {
        for ( j = 0; j < ((*space_size)) ; j++) {
          printf("%lf \t",mat_Dists1[*space_size*i+j]);
        }
	 }
	 printf("\nmat_Dists2\n");
	 for ( i = 0; i < ((*space_size)) ; i++) {
        for ( j = 0; j < ((*space_size)) ; j++) {
          printf("%lf \t",mat_Dists2[*space_size*i+j]);
        } 
	 } */
	
	/*printf("\n\noccupied\n");
	for ( i = 0; i < (*space_size) ; i++) {
		printf("%d \t",(occupied[i]));
	} */
	
	/*printf("\n\nsigma\n");
	for ( i = 0; i < (*space_dim) ; i++) {
		printf("%lf \t",(sigma[i]));
	}
	
	printf("\n\nlambda\n");
	printf("%lf",*lambda);

	printf("\n\ntau\n");
	printf("%lf",*tau);
	
	printf("\n\nmig\n");
	for ( i = 0; i < (*length_mig) ; i++) {
		printf("%lf \t",(mig[i]));
	}
	
	printf("\n\nmig_event\n");
	for ( i = 0; i < (4+(*space_dim)) ; i++) {
		printf("%lf \t",(mig_event[i]));
	} */
	
	// Fonction //
	
	// ##### 1 ##### //
	// initialization departures and destinations //
	long *departures = NULL, *destinations = NULL;
	long length_departures = 0, length_destinations = 0;
	//long length_departures2 = 0;  //-*-
  
	long sum_occupied = 0;
	for (i = 0 ; i < *space_size ; i++) {
		if (occupied[i]>0){
      sum_occupied += 1;
    }
	}
	
	length_departures = sum_occupied ;
	departures = malloc(length_departures * sizeof(long));
	for (i = 0; i < (*space_size); i++) {
		if (occupied[i] > 0 ) {
			//for (j = 0; j < occupied[i]; j++) {
				departures[index] = (i+1);
				index++;	
			//}
		}
	}
	
	// ##### 2 ##### //
	length_destinations = *space_size ; 
  //length_departures2 = *space_size ;//-*-
	destinations = malloc(*space_size * sizeof(long));
  //departures2 = malloc(*space_size * sizeof(long));//-*-
	for (i = 0 ; i < *space_size ; i++) {
  	destinations[i] = (i+1);
		//departures2[i] = (i+1);//-*-
	}
  
	
	/*printf("\n\ndepartures\n");
	for ( i = 0 ; i < length_departures ; i++) {
		printf("[%d] %ld\n", i, departures[i]);
	}
	
	printf("\n\ndestinations\n");
	for ( i = 0 ; i < length_destinations ; i++) {
		printf("[%d] %ld\n", i, destinations[i]);
	} */
	
	// ##### 3 ##### //
	// Distance matrix //
	long dim_mat_dists = length_departures * length_destinations;
  //long dim_mat_dists = length_departures2 * length_destinations;//-*-
/*[for local distance matrix caculation]
  double *mat_distsL = NULL, *mat_distsl = NULL;
	mat_distsL = malloc(dim_mat_dists * sizeof(double));
	mat_distsl = malloc(dim_mat_dists * sizeof(double));
  for (i = 0 ; i < length_departures ; i++) {
	//for (i = 0 ; i < length_departures2 ; i++) {//-*-
		for (j = 0; j < length_destinations ; j++) {
			//mat_distsL[i*length_destinations + j] = fabs(space[departures2[i]-1] - space[destinations[j]-1]);//-*-
			//mat_distsl[i*length_destinations + j] = fabs(space[(*space_size+departures2[i]-1)] - space[(*space_size+destinations[j]-1)]);//-*-
  		//mat_distsL[i*length_destinations + j] = fabs(space[departures[i]-1] - space[destinations[j]-1]);                              // euclidean
			//mat_distsl[i*length_destinations + j] = fabs(space[(*space_size+departures[i]-1)] - space[(*space_size+destinations[j]-1)]);  // euclidean
    //recup la good valeur dans la mtrice 
      mat_distsL[i*length_destinations + j] = distkm( space[departures[i]-1], space[destinations[j]-1], 
        (space[(*space_size+departures[i]-1)]+space[(*space_size+destinations[j]-1)])/2, (space[(*space_size+departures[i]-1)]+space[(*space_size+destinations[j]-1)])/2 );  // km
			  //printf("space : %f",space)
		//	printf("space[departures[i]-1] : %f", space[departures[i]-1]);
		//	printf("space[destinations[j]-1] : %f", space[destinations[j]-1]);
		//	printf("(space[(*space_size+departures[i]-1)]+space[(*space_size+destinations[j]-1)])/2 : %f", (space[(*space_size+departures[i]-1)]+space[(*space_size+destinations[j]-1)])/2);
		//	printf("(space[(*space_size+departures[i]-1)]+space[(*space_size+destinations[j]-1)])/2 : %f", (space[(*space_size+departures[i]-1)]+space[(*space_size+destinations[j]-1)])/2);
			
			mat_distsl[i*length_destinations + j] = distkm( (space[departures[i]-1]+space[destinations[j]-1])/2, (space[departures[i]-1]+space[destinations[j]-1])/2, 
        space[(*space_size+departures[i]-1)], space[(*space_size+destinations[j]-1)] );  // km
		//	printf("(space[departures[i]-1]+space[destinations[j]-1])/2 : %f", (space[departures[i]-1]+space[destinations[j]-1])/2);
		//	printf("(space[departures[i]-1]+space[destinations[j]-1])/2 : %f", (space[departures[i]-1]+space[destinations[j]-1])/2);
		//	printf("space[(*space_size+departures[i]-1)] : %f", space[(*space_size+departures[i]-1)];
		//	printf("space[(*space_size+destinations[j]-1)] : %f", space[(*space_size+destinations[j]-1)];
		}
	}*/
  
  

	 /*
	 printf("\n\nmat_distsL\n");
	 for (i = 0; i < dim_mat_dists; i++) {
	   printf("%lf , ",mat_distsL[i]);
	 }
	 printf("\n\nmat_distsl\n");
	 for (i = 0; i < dim_mat_dists; i++) {
	   printf("%lf , ",mat_distsl[i]);
	 } */
	
	// ##### 4 ##### //
	// Density matrix //
	double mean = 0, sdL, sdl ;
  sdL = sqrt(sigma[0]);
  sdl = sqrt(sigma[1]);
	int b_log = 1;
/*[for local distance matrix caculation]
 double *densitymat = NULL;
	densitymat = malloc(dim_mat_dists * sizeof(double));
printf("\n");
	for (i = 0; i < length_destinations; i++) {
  	for (j = 0; j < length_departures; j++) {
  	  printf("%ld--",i + (length_destinations * j));
		//for (j = 0; j < length_departures2; j++) {//-*-
  		//densitymat[i + (length_destinations * j)] = ((dnorm(mat_distsL[i + (length_destinations * j)] , mean , sdL , b_log)) + (dnorm(mat_distsl[i + (length_destinations * j)] , mean , sdl , b_log))) - ( log((pnorm(1, mean, sdL,1,0)- pnorm(0, mean, sdL,1,0))) + log((pnorm(1,mean,sdl,1,0) - pnorm(0,mean,sdl,1,0)))) ;
			densitymat[i + (length_destinations * j)] = exp( (dnorm(mat_distsL[i + (length_destinations * j)] , mean , sdL , b_log) - dnorm(0 , mean , sdL , b_log)) +
                                                       (dnorm(mat_distsl[i + (length_destinations * j)] , mean , sdl , b_log) - dnorm(0 , mean , sdl , b_log)) )
                                                  / length_destinations ;
		}
  }
printf("\n");*/
	// Alternative: Density matrix using input distance matrices//
/* */	double *densitymat = NULL;
	densitymat = malloc(dim_mat_dists * sizeof(double));
	int k;
	//printf("dim=%ld\n",dim_mat_dists);
	for (i = 0; i < *space_size; i++) {
	  k=0;
	  for (j = 0; j < *space_size; j++) {
      if (occupied[j]>0){
  	    //printf("%d--",k);
  	    densitymat[i + (length_destinations * (j-k))] = exp( (dnorm(mat_Dists1[i + (length_destinations * j)] , mean , sdL , b_log) - dnorm(0 , mean , sdL , b_log)) +
                                                             (dnorm(mat_Dists2[i + (length_destinations * j)] , mean , sdl , b_log) - dnorm(0 , mean , sdl , b_log)) )
                                                        / length_destinations ;
      }else{ // keep track of the unoccupied locations to get correct index
        k++;
      }
    }
	}
	
	 
   //if (sum_occupied==1){
     /*printf("\n\nF_{i,j}\n");
  	 for (i = 0; i < dim_mat_dists ; i++) {
  	   printf("%e , ",densitymat[i]);
  	 }
  	 printf("\nF_{i,j} calculated with input distance matrices\n");
  	 for (i = 0; i < dim_mat_dists ; i++) {
  	   printf("%e , ",densitymati[i]);
  	 }
  	 printf("\n\n"); */
   //}
	
	// ##### 5 ##### //
  
  /* printf("\n\noccupied\n");
	for ( i = 0; i < (*space_size) ; i++) {
		printf("%lf \t",(double)occupied[i]);
	}*/
  
  //FILE *fp;
  //if (sum_occupied==1){fp = fopen("results.dat", "w");}
	// Lambda //
	double *L = NULL;
	L = malloc(length_destinations * sizeof(double));
	for (i = 0 ; i < length_destinations ; i++) {
		L[i] = (double)occupied[destinations[(i)]-1];
		if (L[i]==0) {
			L[i] = 1;
		} else if (L[i]>0) {
			L[i] = *lambda ;
		}
	}
	/*printf("\n\nL\n");
	for (i = 0; i < length_destinations ; i++) {
	  printf("%lf , ",L[i]);
    //printf("%d , ",occupied[i]);
	} */
	
	// ##### 6 ##### //
	// R //
  double *R= NULL;
  R = malloc(dim_mat_dists * sizeof(double));
  double Rs=0;
  for ( i = 0 ; i < length_departures ; i++ ) {
	//for ( i = 0 ; i < length_departures2 ; i++ ) {//-*-
		for ( j =0 ; j < length_destinations ; j++) {
			densitymat[i*length_destinations+j] = densitymat[i*length_destinations+j] * L[j];
      R[i*length_destinations+j] = densitymat[i*length_destinations+j] * *tau;
      //if (occupied[destinations[(j)]-1]>=1) Rs += (R[i*length_destinations+j] * (double)occupied[destinations[(j)]-1]);
      if (occupied[departures[(i)]-1]>=1) Rs += (R[i*length_destinations+j] * (double)occupied[departures[(i)]-1]);
      //printf("%lf ",R[i*length_destinations+j]);
      //if (sum_occupied==1){fprintf(fp,"%d-%d, %lf, %ld, %d\n",i,j,L[j],destinations[(j)],occupied[destinations[(j)]-1]);}
		}
    //printf("\n");
	}
  //if (sum_occupied==1){fclose(fp);}


   /*
   if (sum_occupied==1){
     printf("\n\nR_{i,j}\n");
  	 for (i = 0; i < dim_mat_dists ; i++) {
  	   printf("%e , ",R[i]);
  	 }
   } */
   
	 /*printf("\n\ndensitymat6\n");
	 for (i = 0; i < dim_mat_dists ; i++) {
	   printf("%lf , ",densitymat[i]);
	 } */
	
	// ##### 7 ##### //
  /*double sumd=0;
	for (i = 0 ; i < dim_mat_dists ; i++) {
		if (IsFiniteNumber(densitymat[i])==0) {
			densitymat[i]=0;
		}
    sumd += densitymat[i];
	}*/
  if (Rs==0) {
    /*printf("\nno possible move\n");
    printf("\n\ndensitymat7\n");
    for (i = 0; i < dim_mat_dists ; i++) {
	     printf("%lf , ",R[i]);
	  }*/
  	mig_event[0] = 0;
		mig_event[1] = 0;
		mig_event[2] = 0;
		mig_event[3] = INFINITY;
		for (i = 0; i < *space_dim ; i++) {
			mig_event[4+i] = 0;
		}
		return;
  }

	 /*printf("\n\ndensitymat7\n");
	 for (i = 0; i < dim_mat_dists ; i++) {
	 printf("%lf , ",densitymat[i]);
	 } */
   	
	double RowSums = 0;
	double *p = NULL;
  p = malloc(length_departures*sizeof(double));
	for ( i =0 ; i < length_departures ; i++ ) {
		RowSums = 0;
		for ( j = 0 ; j < length_destinations ; j++) {
			RowSums += R[(length_destinations*i)+j];
		}
		p[i] = RowSums/Rs ;
	}
	/*printf("\n\np\n");
	 for ( i = 0 ; i < length_departures ; i++) {
	 printf("%lf,",p[i]);
	 } */
	 
	/*printf("\n\nexp_param\n");
	for ( i = 0 ; i < length_departures ; i++) {
	printf("%lf,",exp_param[i]);
	} */
	
  /*printf("\n\nsum_exp_param\n");
	printf("%lf,",sum_exp_param); */
	
	double proba_event = 0;
	double wait_time = 0 ;
  int lstart = 0 , lgoto = 0 ;
	
	// ##### 18 ##### //
	if ( *length_mig > 1) { // find the migration event: from mig[1] to mig[2] and return its probability
  	for ( i = 0 ; i < length_departures ; i++ ) {
			if ( departures[i]==(long)mig[0]) {
				lstart = i;
//printf("x0=%f",space[departures[i]-1]);
//printf(",y0=%f ",space[(*space_size+departures[i]-1)]);
				break;
			}
		}
		for ( i = 0 ; i < length_destinations ; i++ ) {
			if ( destinations[i]==(long)mig[1]) {
				lgoto = i;
//printf("x1=%f",space[destinations[i]-1]);
//printf(",y1=%f ,",space[(*space_size+destinations[i]-1)]);
				break;
			}
		}
		// ##### 19 ##### //
    if (sum_occupied==1){ //root
      double sum_row=0;
      for ( i = 0 ; i < length_destinations ; i++ ) {
        sum_row += R[(length_destinations*lstart)+i];
      }
//printf("lstart=%i\n",lstart);
//printf("lgoto=%i\n",lgoto);
//printf("sum_row=%f\n",sum_row);
//printf("R=%f\n",R[(length_destinations*(long)(lstart))+(long)(lgoto)]);
		  proba_event = R[(length_destinations*(long)(lstart))+(long)(lgoto)] / sum_row;
    }else{
      proba_event = R[(length_destinations*(long)(lstart))+(long)(lgoto)] * exp(-Rs*mig[2]);
    }
		wait_time = 0;
		/*printf("\n\ndexp_proba\n");
  	printf("%lf",dexp_proba(mig[2],sum_exp_param)); */
//printf("from %.0f to %.0f, waiting_time=%f, Rs=%e, proba_event=%e (F_{i,j}=%e,R_{i,j}=%e)\n",mig[0],mig[1],mig[2],Rs,proba_event,densitymat[(length_destinations*(long)(lstart))+(long)(lgoto)],R[(length_destinations*(long)(lstart))+(long)(lgoto)]);
	}
	// ##### 21 ##### //
	else { // sample one migration event
    //double sum_densitymat = 0;
    //for (i = 0; i < dim_mat_dists ; i++) {
    //  sum_densitymat += densitymat[i];
	  //}
    wait_time = rexp_proba(Rs);
		/*printf("\n\nwait_time\n");
		printf("%lf",wait_time);*/
    double sum_p = 0;
  	for ( i = 0 ; i < length_departures ; i++) {
  		sum_p += p[i];
  	}
    
    int *x;
    x = malloc(length_departures*sizeof(int));
  	double *probax;
		probax = malloc(length_departures*sizeof(double));
    for ( i = 0 ; i < length_departures ; i++ ) {
    	x[i] = i;
      probax[i] = p[i]/sum_p;
      //printf("%lf\n",probax[i]);
    }
    lstart = sample_proba(x, probax, length_departures);
    //printf("lstart: %d\n",lstart);
    
    int *y;
    y = malloc(length_destinations*sizeof(int));
    double *probay;
		probay = malloc(length_destinations*sizeof(double));
    double sum_row=0;
    for ( i = 0 ; i < length_destinations ; i++ ) {
      sum_row += R[(length_destinations*lstart)+i];
    }
    //printf("sum_row %lf\n",sum_row);
    for ( i = 0 ; i < length_destinations ; i++ ) {
      y[i] = i;
      probay[i] = R[(length_destinations*lstart)+i]/sum_row;
      //printf("[%d] %lf\n",i,probay[i]);
    }
    lgoto = sample_proba(y, probay, length_destinations);
    /*int *x, event;
  	x = malloc(dim_mat_dists*sizeof(int));
		double *proba;
		proba = malloc(dim_mat_dists * sizeof(double));
    for ( i = 0 ; i < length_departures ; i++ ) {
    	for ( j =0 ; j < length_destinations ; j++) {
  		  x[(length_destinations*i)+j] = (length_destinations*i)+j;
			  proba[(length_destinations*i)+j] = densitymat[(length_destinations*i)+j]/sum_densitymat ; // need proper proba (sum to 1) to call sample_proba()
  		}
  	}*/
		//////////// SAMPLE ///////////////////
  /*printf("\n\nproba\n");
	for ( i = 0 ; i < dim_mat_dists ; i++) {
	printf("%lf,",proba[i]);
	} */
		//event = sample_proba(x, proba, dim_mat_dists);
    //lstart = event/length_destinations;
    //lgoto = event%length_destinations;
    //proba_event = probax[lstart] * probay[lgoto];
  	proba_event = 0 ; //DEBUGGING return the likelihood of the migration event:
    if (sum_occupied==1){ //root
  	  proba_event = R[(length_destinations*(long)(lstart))+(long)(lgoto)] / sum_row;
    }else{
      proba_event = R[(length_destinations*(long)(lstart))+(long)(lgoto)] * exp(-Rs*wait_time);
    }
		free(x);
  	free(y);
		free(probax);
		free(probay);
	}
	/*printf("\n\nlstart \n %d",lstart);
	 printf("\n\nlgoto \n %d",lgoto);
	 printf("\n\nproba_event \n %lf", proba_event);
	 printf("\n\nwait_time \n %lf", wait_time); */
	
	
	// ##### 22 ##### //
	mig_event[0] = departures[(long)(lstart)]  ;
	mig_event[1] = destinations[(long)(lgoto)]  ;
	mig_event[2] = proba_event ;
	mig_event[3] = wait_time ;
	/* this is wrong if distkm:
	for ( i = 0 ; i < *space_dim ; i++) {
		mig_event[4+i] = fabs(space[((*space_size)*i+(destinations[(long)(lgoto)]-1))] - space[((*space_size)*i+(departures[(long)(lstart)]-1))]);
	}*/
	mig_event[4] = mat_Dists1[(destinations[(long)(lgoto)]-1) + (length_destinations * (departures[(long)(lstart)]-1))] ;
	mig_event[5] = mat_Dists2[(destinations[(long)(lgoto)]-1) + (length_destinations * (departures[(long)(lstart)]-1))] ;

	
	/*printf("\n\nmig_event\n");	
	for ( i = 0 ; i < 6 ; i++) {
	 printf("%lf  ",mig_event[i]);
	} */
	
	free(departures);
	free(destinations);
/*[for local distance matrix caculation]	
	free(mat_distsL);
	free(mat_distsl);
*/
	free(densitymat);
	free(L);
	free(R);
  free(p);
	
	return;
	
}


int IsFiniteNumber(double x) 
{
	return (x <= DBL_MAX && x >= -DBL_MAX); 
}


long sample_proba( int x[], double proba[], long length)
{
	int i = 0, res = 0;
  double a, sum_proba = 0;
	//double a = (random()/(double)RAND_MAX);
  GetRNGstate();
  a = unif_rand();
  PutRNGstate();
	for ( i = 0 ; i < length ; i++) {
		sum_proba += proba[i];
		if (a < sum_proba) {
			res = x[i];
			break;
		}
	}
	return res;
}

double rexp_proba(double lambda_exp)
{
	double res = 0 ;
	//res = -log((random()/(double)RAND_MAX))/(lambda_exp);
  GetRNGstate();
  res = -log(unif_rand())/(lambda_exp);
  PutRNGstate();
	return res;
}


double dexp_proba(double x, double lambda_exp)
{
	double res = 0;
	res = lambda_exp * exp((-lambda_exp) * x);
	return res;
}

double distkm(double lat1, double lat2, double long1, double long2)
{
  double lat1r, lat2r, long1r, long2r, d, R,x,y;
  lat1r = (lat1/180) * M_PI ;
  lat2r = (lat2/180) * M_PI ;
  long1r = (long1/180) * M_PI ;
  long2r = (long2/180) * M_PI ;
  //d=sqrt(((long1-long2)*(long1-long2)) + ((lat1-lat2)*(lat1-lat2)));
  R=6371;
  x=(long2r-long1r)*cos(0.5*(lat2r+lat1r));
  y=lat2r-lat1r;
  d=R*sqrt(x*x+y*y);
  if(d<=100)
  {
    return d;
  }
  else
  {
    if(lat1==lat2 && long1==long2 ){
      return 0;
    }
 
    if((long2r != long1r) && (fabs(long2r-long1r)<0.000001)){
      d = acos( sin(lat1r)*sin(lat2r) + cos(lat1r)*cos(lat2r)) * R ;
      //d=sqrt(((long1-long2)*(long1-long2)) + ((lat1-lat2)*(lat1-lat2)));
    } else{
      d = acos( sin(lat1r)*sin(lat2r) + cos(lat1r)*cos(lat2r)*cos(long2r-long1r)) * R ;
      //d=sqrt(((long1-long2)*(long1-long2)) + ((lat1-lat2)*(lat1-lat2)));
    }
    
    //printf("sin(lat1r):%f sin(lat2r)%f cos(lat1r)%f cos(lat2r):%f produit:%f acos:%f\n",sin(lat1r),sin(lat2r),cos(lat1r),cos(lat2r),(sin(lat1r)*sin(lat2r) + cos(lat1r)*cos(lat2r)), acos( sin(lat1r)*sin(lat2r) + cos(lat1r)*cos(lat2r)));
    //printf("lat1:%f lat2%f long1:%f long2:%f\n",lat1,lat2,long1,long2);
    //printf("lat1r:%f lat2r:%f long1r:%f long2r:%f\n",lat1r,lat2r,long1r,long2r);
    // printf("Dist : %f\n\n",d);
  }
  return d;
}
