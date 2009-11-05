/* $Id$ */
/******************************************************************/
/* This is the wrapper function for to RunCode.R. Its inputs are  */
/* R = number of replicates, P = number of genes, T = number of   */
/* times, K = number of hidden states (0 means no x's estimated), */
/* x0 = initial values of x, yorig = gene expression, a0-sigma0=  */
/* initial values of alpha-sigma, conv1-conv3 = convergence       */
/* criterai, DEst and DvarEst will return the value of the        */
/* posterior mean and variance of D.                              */
/******************************************************************/

#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

/* Load global subprograms */
void MatrixInv(double**, int, double**, double*);
void MatrixMult(double**, int, int, double**, int, double**);
void MatrixSum(double**, double**, double**, int*, int*);
void MatrixTrans(double**, double**, int*, int*);

void RunWrap(int *R, int *P, int *T, int *K, double *xx,
    double *yy, double *alpha, double *beta, double *gamma,
    double *delta, double *v, double *mu, double *sigma,
    double *conv1, double *conv2, double *conv3,
    double *APost, double *BPost, double *CPost, double *DPost,
    double *CvarPost, double *DvarPost, int *alliterations, int *maxiterations,
    int *subiterations)
{

    /* Load subprograms */
    void RunCode(int*, int*, int*, int*, double*, double*, double*,
        double*, double*, double*, double*, double*, double*,
        double*, double*, double*, double*, double*, double*, double*, double*, double*, int*, int*, int*);

    RunCode(R, P, T, K, xx, yy, alpha, beta, gamma, delta, v, mu,
        sigma, conv1, conv2, conv3, APost, BPost, CPost, DPost, CvarPost, DvarPost, alliterations, maxiterations,
        subiterations);

}

/* Include other programs */

/******************************************************************/
/*  This is the counterpart of Overall_Posterior_Mean.R. It has   */
/*  inputs alpha, beta, gamma, delta, v, x, y, K, P, T, and R     */
/*  and returns the value of the overall posterior mean of A, B,  */
/*  C, and D (in the case where x's are estimated), or the        */
/*  overall posterior mean and variance of D (in the case where   */
/*  x's are not estimated).                                       */
/******************************************************************/

void PostMeanOverall(double *alpha, double *beta, double *gamma, double *delta,
	double *v, double ***x, double ***y, int *K, int *P, int *T,
	int *R, double *AA, double *BB, double *CC, double *DD, double *CCvar, double *DDvar)
{

	int i, ii, iii, *all, j, jj, index, *KK, index2;
    double **MNLNinv, **LNMNinv, **JGFGinv, **FGJGinv, ***HNLS,
            ***SNMH, ***EGFQ, ***QGJE, **matk1, **matp1, **Dvartemp,
            **D, ***Dvar, ***Cvar, **C, **B, **A;

	/* Load subprograms */

	void SimplifyNoX(double*, double*, double***, int*, int*, int*, int*,
        	double**, double**);
	void SimplifyX(double*, double*, double*, double*, double*, double***,
		double***, int*, int*, int*, int*, int*, double**, double**, double**,
        	double**, double***, double***, double***,double***);

    /* Allocate memory */

	all = (int*) calloc(1, sizeof(int));
	KK = (int*) calloc(1, sizeof(int));
	*all = 1;
	*KK = 1;
	if(*K > 0)
	{
	    *KK = *K;
	}

	A = (double**) calloc(*KK, sizeof(double*));
	B = (double**) calloc(*KK, sizeof(double*));
    	C = (double**) calloc(*P, sizeof(double*));
	Cvar = (double***) calloc(*P, sizeof(double**));
	D = (double**) calloc(*P, sizeof(double*));
	Dvar = (double***) calloc(*P, sizeof(double**));
    for(j=0; j<*KK; j++)
	{
	    *(A+j) = (double*) calloc(*KK, sizeof(double));
	    *(B+j) = (double*) calloc(*P, sizeof(double));
	}
	for(i=0; i<*P; i++)
	{
        *(C+i) = (double*) calloc(*KK, sizeof(double));
	    *(D+i) = (double*) calloc(*P, sizeof(double));
	    *(Dvar+i) = (double**) calloc(*P, sizeof(double*));
	    *(Cvar+i) = (double**) calloc(*KK, sizeof(double*));
	    for(ii=0; ii<*P; ii++)
	    {
	        *(*(Dvar+i)+ii) = (double*) calloc(*P, sizeof(double));
	    }
	   for(j=0; j<*KK; j++)
	   {
		*(*(Cvar+i)+j) = (double*) calloc(*KK, sizeof(double));
	   }
	}

	/* Begin function -- no x's */

    Dvartemp = (double**) calloc(*P, sizeof(double*));
    for(i=0; i<*P; i++)
    {
        *(Dvartemp+i) = (double*) calloc(*P, sizeof(double));
    }

	if(*K==0)
	{
        SimplifyNoX(delta, v, y, P, T, R, all, D, Dvartemp);
        for(i=0; i<*P; i++)
        {
            for(ii=0; ii<*P; ii++)
            {
                for(iii=0; iii<*P; iii++)
                {
                    *(*(*(Dvar+i)+ii)+iii) = (1.0/(*(v+i))) * (*(*(Dvartemp+ii)+iii));
                }
            }
        }
	}

    /* Free memory */
    for(i=0; i<*P; i++)
    {
        free(*(Dvartemp+i));
    }
    free(Dvartemp);

	/* Begin function -- x's */
    /* Allocate MNLNinf, LNMNinf, JGFGinv, FGJGinv, HNLS, SNMH, EGFQ, QFJE */
    MNLNinv = (double**) calloc(*KK, sizeof(double*));
    LNMNinv = (double**) calloc(*P, sizeof(double*));
    JGFGinv = (double**) calloc(*KK, sizeof(double*));
    FGJGinv = (double**) calloc(*P, sizeof(double*));
    HNLS = (double***) calloc(*KK, sizeof(double**));
    SNMH = (double***) calloc(*KK, sizeof(double**));
    EGFQ = (double***) calloc(*P, sizeof(double**));
    QGJE = (double***) calloc(*P, sizeof(double**));
    for(j=0; j<*KK; j++)
    {
        *(MNLNinv+j) = (double*) calloc(*KK, sizeof(double));
        *(JGFGinv+j) = (double*) calloc(*KK, sizeof(double));
        *(HNLS+j) = (double**) calloc(*KK, sizeof(double*));
        *(SNMH+j) = (double**) calloc(*P, sizeof(double*));
        for(jj=0; jj<*KK; jj++)
        {
            *(*(HNLS+j)+jj) = (double*) calloc(1, sizeof(double));
        }
        for(i=0; i<*P; i++)
        {
            *(*(SNMH+j)+i) = (double*) calloc(1, sizeof(double));
        }
    }
    for(i=0; i<*P; i++)
    {
        *(LNMNinv+i) = (double*) calloc(*P, sizeof(double));
        *(FGJGinv+i) = (double*) calloc(*P, sizeof(double));
        *(EGFQ+i) = (double**) calloc(*KK, sizeof(double*));
        *(QGJE+i) = (double**) calloc(*P, sizeof(double*));

        for(j=0; j<*KK; j++)
        {
            *(*(EGFQ+i)+j) = (double*) calloc(1, sizeof(double));
        }
        for(ii=0; ii<*P; ii++)
        {
            *(*(QGJE+i)+ii) = (double*) calloc(1, sizeof(double));
        }
    }
    matk1 = (double**) calloc(*KK, sizeof(double*));
    matp1 = (double**) calloc(*P, sizeof(double*));
    for(j=0; j<*KK; j++)
    {
        *(matk1+j) = (double*) calloc(1, sizeof(double));
    }
    for(i=0; i<*P; i++)
    {
        *(matp1+i) = (double*) calloc(1, sizeof(double));
    }

	if(*K > 0)
	{
        SimplifyX(alpha, beta, gamma, delta, v, x, y, K, P, T, all, R,
            MNLNinv, LNMNinv, JGFGinv, FGJGinv, HNLS, SNMH, EGFQ,QGJE);

        for(j=0; j<*K; j++)
        {
            MatrixMult(MNLNinv, *K, *K, *(HNLS+j), 1, matk1);
            MatrixMult(LNMNinv, *P, *P, *(SNMH+j), 1, matp1);
            for(jj=0; jj<*K; jj++)
            {
                *(*(A+j)+jj) = *(*(matk1+jj));
            }
            for(i=0; i<*P; i++)
            {
                *(*(B+j)+i) = *(*(matp1+i));
            }
        }

        for(i=0; i<*P; i++)
        {
            MatrixMult(FGJGinv, *P, *P, *(QGJE+i), 1, matp1);
            MatrixMult(JGFGinv, *K, *K, *(EGFQ+i), 1, matk1);
            for(ii=0; ii<*P; ii++)
            {
                *(*(D+i)+ii) = *(*(matp1+ii));
                for(iii=0; iii<*P; iii++)
                {
                    *(*(*(Dvar+i)+ii)+iii) = (1.0/(*(v+i))) * (*(*(FGJGinv+ii)+iii));
                }
            }
            for(j=0; j<*K; j++)
            {
               *(*(C+i)+j) =  *(*(matk1+j));
		for(jj=0; jj<*K; jj++)
		{
			*(*(*(Cvar+i)+j)+jj) = (1.0/(*(v+i))) * (*(*(JGFGinv+j)+jj));
		}
            }
        }
	}

	/* Set DD = D, CC = C, CCvar = Cvar, and DDvar = Dvar */
	index = 0;
	index2 = 0;
	for(i=0; i<*P; i++)
	{
	    for(ii=0; ii<*P; ii++)
	    {
	        for(iii=0; iii<*P; iii++)
	        {
	           *(DDvar+index) = *(*(*(Dvar+i)+ii)+iii);
	           index++;
	        }
	    }
	    for(j=0; j<*K; j++)
	    {
		for(jj=0; jj<*K; jj++)
		{
			*(CCvar+index2) = *(*(*(Cvar+i)+j)+jj);
			index2++;
		}
	    }
	}

	index = 0;
	index2 = 0;
	for(i=0; i<*P; i++)
	{
	    for(ii=0; ii<*P; ii++)
	    {
	        *(DD + index) = *(*(D+i)+ii);
	        index++;
	    }
	    for(j=0; j<*K; j++)
	    {
	        *(CC + index2) = *(*(C+i)+j);
		index2++;
	    }
	}
	index = 0;
    index2 = 0;
    for(j=0; j<*K; j++)
    {
        for(jj=0; jj<*K; jj++)
        {
            *(AA + index) = *(*(A+j)+jj);
            index++;
        }
        for(i=0; i<*P; i++)
        {
            *(BB + index2) = *(*(B+j)+i);
            index2++;
        }
    }

	/* Free memory */
	/* If x's */
    for(j=0; j<*KK; j++)
    {
        for(jj=0; jj<*KK; jj++)
        {
            free(*(*(HNLS+j)+jj));
        }
        for(i=0; i<*P; i++)
        {
            free(*(*(SNMH+j)+i));
        }
        free(*(HNLS+j));
        free(*(SNMH+j));
        free(*(MNLNinv+j));
        free(*(JGFGinv+j));
    }

    for(i=0; i<*P; i++)
    {
        for(j=0; j<*KK; j++)
        {
            free(*(*(EGFQ+i)+j));
        }
        for(ii=0; ii<*P; ii++)
        {
            free(*(*(QGJE+i)+ii));
        }
        free(*(EGFQ+i));
        free(*(QGJE+i));
        free(*(LNMNinv+i));
        free(*(FGJGinv+i));
    }
    free(HNLS);
    free(SNMH);
    free(MNLNinv);
    free(JGFGinv);
    free(EGFQ);
    free(QGJE);
    free(LNMNinv);
    free(FGJGinv);

    for(j=0; j<*KK; j++)
    {
        free(*(matk1+j));
    }
    for(i=0; i<*P; i++)
    {
        free(*(matp1+i));
    }
    free(matk1);
    free(matp1);

	/* Release memory of D's */
    for(j=0; j<*KK; j++)
	{
	    free(*(A+j));
	    free(*(B+j));
	}
	for(i=0; i<*P; i++)
	{
	    for(ii=0; ii<*P; ii++)
	    {
	        free(*(*(Dvar+i)+ii));
	    }
	    for(j=0; j<*KK; j++)
	    {
		free(*(*(Cvar+i)+j));
	    }
	    free(*(Cvar+i));
	    free(*(Dvar+i));
	    free(*(D+i));
        free(*(C+i));
	}
	free(Cvar);
	free(Dvar);
	free(D);
    	free(C);
	free(KK);
	free(all);
    	free(A);
    	free(B);
}




/******************************************************************/
/* This is counterpart to RunCode.R.  It takes as inputs          */
/* R = number of replicates, P = number of genes, T = number of   */
/* times, K = number of hidden states (0 means no x's estimated), */
/* x0 = initial values of x, yorig = gene expression, a0-sigma0=  */
/* initial values of alpha-sigma, conv1-conv3 = convergence       */
/* criteria, DPost and DvarPost will return the value of the      */
/* posterior mean and variance of D.                              */
/* ****************************************************************/

void RunCode(int *R, int *P, int *T, int *K, double *xx,
    double *yy, double *alpha, double *beta, double *gamma,
    double *delta, double *v, double *mu, double *sigma,
    double *conv1, double *conv2, double *conv3,
    double *APost, double *BPost, double *CPost, double *DPost, double *CvarPost,
    double *DvarPost, int *alliterations, int *maxiterations, int *subiterations)
{
    int i, j, r, t, index, *KK;
    double ***x, ***y;

    /* Load subprograms */
    void FullAlgorithm(double***, double***, double*, double*,
        double*, double*, double*, double*, double*, int*, int*,
        int*, int*, double*, double*, double*, int*, int*, int*);
    void PostMeanOverall(double*, double*, double*, double*,
        double*, double***, double***, int*, int*, int*, int*, double*, double*,
        double*, double*, double*, double*);

    /* Allocate y, read in data */
    y = (double***) calloc(*R, sizeof(double**));
    KK = (int*) calloc(1, sizeof(int));
    index = 0;
    *KK = 1;
    if(*K > 0) {
        *KK = *K;
    }
    for(r=0; r<*R; r++)
    {
        *(y+r) = (double**) calloc(*P, sizeof(double*));
        for(i=0; i<*P; i++)
        {
            *(*(y+r)+i) = (double*) calloc(*T, sizeof(double));
            for(t=0; t<*T; t++)
            {
                *(*(*(y+r)+i)+t) = *(yy+index);
                index++;
            }
        }
    }

    /* Allocate x's if K > 0, read in data */
    x = (double***) calloc(*R, sizeof(double**));
    for(r=0; r<*R; r++)
    {
        *(x+r) = (double**) calloc(*KK, sizeof(double*));
        for(j=0; j<*KK; j++)
        {
            *(*(x+r)+j) = (double*) calloc(*T, sizeof(double));
        }
    }

    index = 0;
    if(*K > 0)
    {
        for(r=0; r<*R; r++)
        {
            for(j=0; j<*K; j++)
            {
                for(t=0; t<*T; t++)
                {
                    *(*(*(x+r)+j)+t) = *(xx+index);
                    index++;
                }
            }
        }
    }


    /*****************************/
    /*  Run Full algorithm       */
    /*****************************/

    FullAlgorithm(y, x, alpha, beta, gamma, delta, v, mu, sigma,
        K, P, R, T, conv1, conv2, conv3, alliterations, maxiterations, subiterations);
   /* Rprintf("\n"); */

    /****************************************/
    /*  Find posterior mean & variance      */
    /****************************************/

    PostMeanOverall(alpha, beta, gamma, delta, v, x, y, K, P, T, R, APost, BPost, CPost, DPost, CvarPost, DvarPost);

    /* Read in final estimate of x's */
    if(*K > 0)
    {
        index = 0;
        for(r=0; r<*R; r++)
        {
            for(j=0; j<*KK; j++)
            {
                for(t=0; t<*T; t++)
                {
                    *(xx+index) = *(*(*(x+r)+j)+t);
                    index++;
                }
            }
        }
    }

    /* Release memory of x's if K > 0 */
    for(r=0; r<*R; r++)
    {
        for(j=0; j<*KK; j++)
        {
            free(*(*(x+r)+j));
        }
        free(*(x+r));
    }
    free(x);

    /* Release memory of y's */
    for(r=0; r<*R; r++)
    {
        for(i=0; i<*P; i++)
        {
            free(*(*(y+r)+i));
        }
        free(*(y+r));
    }
    free(y);
    free(KK);
}



void FullAlgorithm(double ***y, double ***x, double *alpha,
    double *beta, double *gamma, double *delta, double *v,
    double *mu, double *sigma, int *K, int *P, int *R, int *T,
    double *conv1, double *conv2, double *conv3, int *alliterations, int *maxiterations,
    int *subiterations)
{
    int iter=1, j, i, r, *KK;
    double overalldiff, ***A, ***B, ***C, ***D, ***Dvar,
        *alpha0, *beta0, *gamma0, *delta0, *v0,
        *alphaold, *betaold, *gammaold, *deltaold, *vold;
    double alphadiff, betadiff, gammadiff, deltadiff, sumnum, sumden;

    /* Load subprograms */
    void EmTypeConv(double*, double*, double*, double*, double*,
        double***, double***, int*, int*, int*, int*, double*, double*, int*);
    void PostMeanR(double*, double*, double*, double*, double*, double***,
        double***, int*, int*, int*, int*, double***, double***,
        double***, double***, double***);
    void Kalman(double***, double***, double***, double***, double***,
        double*, double*, double*, int*, int*, int*, int*, double***);


    /* Allocate memory */
    KK = (int*) calloc(1, sizeof(int));
    *KK = 1;
    if(*K > 0) {
        *KK = *K;
    }

    A = (double***) calloc(*R, sizeof(double**));
    B = (double***) calloc(*R, sizeof(double**));
    C = (double***) calloc(*R, sizeof(double**));
    D = (double***) calloc(*R, sizeof(double**));
    Dvar = (double***) calloc(*R, sizeof(double**));
    for(r=0; r<*R; r++)
    {
        *(A+r) = (double**) calloc(*KK, sizeof(double*));
        *(B+r) = (double**) calloc(*KK, sizeof(double*));
        *(C+r) = (double**) calloc(*P, sizeof(double*));
        *(D+r) = (double**) calloc(*P, sizeof(double*));
        *(Dvar+r) = (double**) calloc(*P, sizeof(double*));
        for(j=0; j<*KK; j++)
        {
            *(*(A+r)+j) = (double*) calloc(*KK, sizeof(double));
            *(*(B+r)+j) = (double*) calloc(*P, sizeof(double));
        }
        for(i=0; i<*P; i++)
        {
            *(*(C+r)+i) = (double*) calloc(*KK, sizeof(double));
            *(*(D+r)+i) = (double*) calloc(*P, sizeof(double));
            *(*(Dvar+r)+i) = (double*) calloc(*P, sizeof(double));
        }
    }
    alpha0 = (double*) calloc(*KK, sizeof(double));
    beta0 = (double*) calloc(*P, sizeof(double));
    gamma0 = (double*) calloc(*KK, sizeof(double));
    delta0 = (double*) calloc(*P, sizeof(double));
    v0 = (double*) calloc(*P, sizeof(double));
    alphaold = (double*) calloc(*KK, sizeof(double));
    betaold = (double*) calloc(*P, sizeof(double));
    gammaold = (double*) calloc(*KK, sizeof(double));
    deltaold = (double*) calloc(*P, sizeof(double));
    vold = (double*) calloc(*P, sizeof(double));


    /* If no x's */
    if(*K == 0)
    {
        EmTypeConv(alpha, beta, gamma, delta, v, x, y, K, P, T, R, conv1, conv2, subiterations);
    }

    /* If x estimates */
    if(*K > 0)
    {
        /* Three convergence criteria:
           Convergence for first sub-loop of EM (pre-v) is conv1
           Convergence for second sub-loop of EM (post-v) is conv2
           Overall algorithm convergence is conv3 */
        /*********************************************************/
        /*  Initial hyperparameter estimates, based on initial x */
        /*  and initial hyperparameters                          */
        /*  New values are stored back in alpha, beta, etc.      */
        /*********************************************************/
        for(j=0; j<*K; j++)
        {
            *(alpha0+j) = *(alpha+j);
            *(gamma0+j) = *(gamma+j);
        }
        for(i=0; i<*P; i++)
        {
            *(beta0+i) = *(beta+i);
            *(delta0+i) = *(delta+i);
            *(v0+i) = *(v+i);
        }
        EmTypeConv(alpha, beta, gamma, delta, v, x, y, K, P, T, R,
            conv1, conv2, subiterations);

        for(j=0; j<*K; j++)
        {
            *(alphaold+j) = *(alpha+j);
            *(gammaold+j) = *(gamma+j);
        }
        for(i=0; i<*P; i++)
        {
            *(betaold+i) = *(beta+i);
            *(deltaold+i) = *(delta+i);
            *(vold+i) = *(v+i);
        }
        Rprintf("0 ");
        overalldiff = 100.0;

        /**************************************************************************************/
        /* Stop if we reach convergence criterion OR if we get up to 100 overall iterations   */
        /* If we are over 100 iterations, we will start over with a new set of initial values */
        /**************************************************************************************/

        while(overalldiff > *conv3)
        {
            if(*alliterations > *maxiterations) break;
            PostMeanR(alpha, beta, gamma, delta, v, x, y, K, P, T,
                R, A, B, C, D, Dvar);

            /******************************************************/
            /*  Kalman filter and smoother estimates of x         */
            /******************************************************/
            Kalman(y, A, B, C, D, v, mu, sigma, K, P, T, R, x);

            /*  Normally we would re-estimate mu and sigma here */

            /******************************************************/
            /*  Running EM algorithm based on new x values        */
            /*  Start with original initial values                */
            /******************************************************/

            for(j=0; j<*K; j++)
            {
                *(alpha+j) = *(alpha0+j);
                *(gamma+j) = *(gamma0+j);
            }
            for(i=0; i<*P; i++)
            {
                *(beta+i) = *(beta0+i);
                *(delta+i) = *(delta0+i);
                *(v+i) = *(v0+i);
            }
            EmTypeConv(alpha, beta, gamma, delta, v, x, y, K, P, T, R,
                conv1, conv2, subiterations);

            /******************************************************/
            /*  Check convergence of hyperparameters              */
            /******************************************************/
            sumnum = 0;
            sumden = 0;
            for(j=0; j<*K; j++)
            {
                sumden += (*(alphaold+j)) * (*(alphaold+j));
            }
            for(j=0; j<*K; j++)
            {
                sumnum += ((*(alpha+j) - *(alphaold+j))*(*(alpha+j) - *(alphaold+j)))/sumden;
            }
            alphadiff = sqrt(sumnum);                                       /* alpha.diff */

            sumnum = 0;
            sumden = 0;
            for(i=0; i<*P; i++)
            {
                sumden += (*(betaold+i)) * (*(betaold+i));
            }
            for(i=0; i<*P; i++)
            {
                sumnum += ((*(beta+i) - *(betaold+i))*(*(beta+i) - *(betaold+i)))/sumden;
            }
            betadiff = sqrt(sumnum);                                       /* beta.diff */

            sumnum = 0;
            sumden = 0;
            for(j=0; j<*K; j++)
            {
                sumden += (*(gammaold+j)) * (*(gammaold+j));
            }
            for(j=0; j<*K; j++)
            {
                sumnum += ((*(gamma+j) - *(gammaold+j))*(*(gamma+j) - *(gammaold+j)))/sumden;
            }
            gammadiff = sqrt(sumnum);                                       /* gamma.diff */

            sumnum = 0;
            sumden = 0;
            for(i=0; i<*P; i++)
            {
                sumden += (*(deltaold+i)) * (*(deltaold+i));
            }
            for(i=0; i<*P; i++)
            {
                sumnum += ((*(delta+i) - *(deltaold+i))*(*(delta+i) - *(deltaold+i)))/sumden;
            }
            deltadiff = sqrt(sumnum);                                       /* delta.diff */

            /* Find maximum convergence criterion */
            overalldiff = alphadiff;
            if(betadiff>overalldiff) {overalldiff = betadiff;}
            if(gammadiff>overalldiff) {overalldiff = gammadiff;}
            if(deltadiff>overalldiff) {overalldiff = deltadiff;}

            /* Set new values = to old values of hyperparameters */
            for(j=0; j<*K; j++)
            {
                *(alphaold+j) = *(alpha+j);
                *(gammaold+j) = *(gamma+j);
            }
            for(i=0; i<*P; i++)
            {
                *(betaold+i) = *(beta+i);
                *(deltaold+i) = *(delta+i);
            }

            /* Update the number of iterations, and loop */
           Rprintf("* %d ", iter);
           *alliterations = iter;
            iter++;
        }

  /*      Rprintf("Estimation complete \n"); */
  }
    /* Free memory */
    for(r=0; r<*R; r++)
    {
        for(j=0; j<*KK; j++)
        {
            free(*(*(A+r)+j));
            free(*(*(B+r)+j));
        }
        for(i=0; i<*P; i++)
        {
            free(*(*(C+r)+i));
            free(*(*(D+r)+i));
            free(*(*(Dvar+r)+i));
        }
        free(*(A+r));
        free(*(B+r));
        free(*(C+r));
        free(*(D+r));
        free(*(Dvar+r));
    }
    free(A);
    free(B);
    free(C);
    free(D);
    free(Dvar);
    free(alpha0);
    free(beta0);
    free(gamma0);
    free(delta0);
    free(v0);
    free(alphaold);
    free(betaold);
    free(gammaold);
    free(deltaold);
    free(vold);
    free(KK);
}



/******************************************************************/
/*  This is the counterpart of Kalman_Filter.R.  It takes as      */
/*  inputs y, the current values of A, B, C, D, v, mu, and sigma  */
/*  (where sigma is a matrix), K, P, T, and R.  It returns the    */
/*  new values of x into the original triple pointer x.           */
/******************************************************************/


void Kalman(double ***y, double ***A, double ***B,
    double ***C, double ***D, double *v, double *mu, double *sigma,
    int *K, int *P, int *T, int *R, double ***x)
{
    int r, t, i, j, ii, jj, index;
    double **yr, **Ar, **Br, **Cr, **Dr, **sigmamat, **xminus,
        **filter, **Pminus, **Pk, **smoother, **Ps;

    /* Load subprograms */
    void KalmanFilter(double**, double**, double**, double**,
        double**, double*, double*, double**, int*, int*, int*,
        double**, double**, double**, double**);

    void KalmanSmoother(double**, double**, double**, double**,
        double**, int*, int*, double**, double**);

    /* Allocate memory */
    yr = (double**) calloc(*P, sizeof(double*));
    Ar = (double**) calloc(*K, sizeof(double*));
    Br = (double**) calloc(*K, sizeof(double*));
    Cr = (double**) calloc(*P, sizeof(double*));
    Dr = (double**) calloc(*P, sizeof(double*));
    sigmamat = (double**) calloc(*K, sizeof(double*));
    xminus = (double**) calloc(*K, sizeof(double*));
    filter = (double**) calloc(*K, sizeof(double*));
    Pminus = (double**) calloc(*K, sizeof(double*));
    Pk = (double**) calloc(*K, sizeof(double*));
    smoother = (double**) calloc(*K, sizeof(double*));
    Ps = (double**) calloc(*K, sizeof(double*));
    for(j=0; j<*K; j++)
    {
        *(Ar+j) = (double*) calloc(*K, sizeof(double));
        *(Br+j) = (double*) calloc(*P, sizeof(double));
        *(sigmamat+j) = (double*) calloc(*K, sizeof(double));
        *(xminus+j) = (double*) calloc(*T, sizeof(double));
        *(filter+j) = (double*) calloc(*T, sizeof(double));
        *(Pminus+j) = (double*) calloc(*K, sizeof(double));
        *(Pk+j) = (double*) calloc(*K, sizeof(double));
        *(smoother+j) = (double*) calloc(*T, sizeof(double));
        *(Ps+j) = (double*) calloc(*K, sizeof(double));
    }
    for(i=0; i<*P; i++)
    {
        *(yr+i) = (double*) calloc(*T, sizeof(double));
        *(Cr+i) = (double*) calloc(*K, sizeof(double));
        *(Dr+i) = (double*) calloc(*P, sizeof(double));
    }

    /* Change sigma into matrix form */
    index = 0;
    for(j=0; j<*K; j++)
    {
        for(jj=0; jj<*K; jj++)
        {
            *(*(sigmamat+j)+jj) = *(sigma+index);
            index++;
        }
    }

    /* Do Kalman filter for each replicate r independently */
    /* Set xminus, filter, Pminus, Pk, smoother, and PS equal to zero at start */

    for(r=0; r<*R; r++)
    {
        /* Read in correct values for y, A, B, C, and D */
        for(i=0; i<*P; i++)
        {
            for(t=0; t<*T; t++)
            {
                *(*(yr+i)+t) = *(*(*(y+r)+i)+t);
            }
            for(j=0; j<*K; j++)
            {
                *(*(Cr+i)+j) = *(*(*(C+r)+i)+j);
            }
            for(ii=0; ii<*P; ii++)
            {
                *(*(Dr+i)+ii) = *(*(*(D+r)+i)+ii);
            }
        }
        for(j=0; j<*K; j++)
        {
            for(jj=0; jj<*K; jj++)
            {
                *(*(Ar+j)+jj) = *(*(*(A+r)+j)+jj);
            }
            for(i=0; i<*P; i++)
            {
                *(*(Br+j)+i) = *(*(*(B+r)+j)+i);
            }
        }

        KalmanFilter(yr, Ar, Br, Cr, Dr, v, mu, sigmamat, K, P,                 /* Kalman filter */
            T, xminus, filter, Pminus, Pk);
        KalmanSmoother(Ar, xminus, filter, Pminus, Pk, K, T,                    /* Kalman smoother */
            smoother, Ps);

        for(j=0; j<*K; j++)
        {
            for(t=0; t<*T; t++)
            {
                *(*(*(x+r)+j)+t) = *(*(smoother+j)+t);                         /* Set x's to smoothed value */
            }
        }
    }

     /* Release memory */
    for(j=0; j<*K; j++)
    {
        free(*(Ar+j));
        free(*(Br+j));
        free(*(sigmamat+j));
        free(*(xminus+j));
        free(*(filter+j));
        free(*(Pminus+j));
        free(*(Pk+j));
        free(*(smoother+j));
        free(*(Ps+j));
    }
    for(i=0; i<*P; i++)
    {
        free(*(yr+i));
        free(*(Cr+i));
        free(*(Dr+i));
    }
    free(yr);
    free(Ar);
    free(Br);
    free(Cr);
    free(Dr);
    free(sigmamat);
    free(xminus);
    free(filter);
    free(Pminus);
    free(Pk);
    free(smoother);
    free(Ps);

}

/* Kalman Filter function */
void KalmanFilter(double **yr, double **Ar, double **Br, double **Cr,
    double **Dr, double *v, double *mu, double **sigmamat, int *K, int *P,
    int *T, double **xminus, double **filter, double **Pminus, double **Pk)
{
    int t, j, jj, i;
    double **gain, **matp1, **matk1, **matk1b, **matp1b, **matkk,
        **xtemp, **ytemp, **xtemp2, **ytemp2, **Art, **matkk2;

    /* Load subprograms */
    void KalmanGain(double**, double**, double*, int*, int*, double**);

    /* Allocate memory */
    gain = (double**) calloc(*K, sizeof(double*));
    matp1 = (double**) calloc(*P, sizeof(double*));
    matp1b = (double**) calloc(*P, sizeof(double*));
    matk1 = (double**) calloc(*K, sizeof(double*));
    matk1b = (double**) calloc(*K, sizeof(double*));
    matkk = (double**) calloc(*K, sizeof(double*));
    xtemp = (double**) calloc(*K, sizeof(double*));
    ytemp = (double**) calloc(*P, sizeof(double*));
    xtemp2 = (double**) calloc(*K, sizeof(double*));
    ytemp2 = (double**) calloc(*P, sizeof(double*));
    Art = (double**) calloc(*K, sizeof(double*));
    matkk2 = (double**) calloc(*K, sizeof(double*));
    for(j=0; j<*K; j++)
    {
        *(gain+j) = (double*) calloc(*P, sizeof(double));
        *(matk1+j) = (double*) calloc(1, sizeof(double));
        *(matk1b+j) = (double*) calloc(1, sizeof(double));
        *(matkk+j) = (double*) calloc(*K, sizeof(double));
        *(xtemp+j) = (double*) calloc(1, sizeof(double));
        *(xtemp2+j) = (double*) calloc(1, sizeof(double));
        *(Art+j) = (double*) calloc(*K, sizeof(double));
        *(matkk2+j) = (double*) calloc(*K, sizeof(double));
    }
    for(i=0; i<*P; i++)
    {
        *(matp1+i) = (double*) calloc(1, sizeof(double));
        *(matp1b+i) = (double*) calloc(1, sizeof(double));
        *(ytemp+i) = (double*) calloc(1, sizeof(double));
        *(ytemp2+i) = (double*) calloc(1, sizeof(double));
    }

    /* Begin Kalman filter */

    for(t=0; t<*T; t++)
    {
        if(t==0)                                                        /* First time point */
        {
            for(j=0; j<*K; j++)
            {
                *(*(xminus+j)) = *(mu+j);                               /* x.minus */
                for(jj=0; jj<*K; jj++)
                {
                    *(*(Pminus+j)+jj) = *(*(sigmamat+j)+jj);               /* P.minus */
                }
            }

            /* Set gain = 0 */
            for(j=0; j<*K; j++)
            {
                for(i=0; i<*P; i++)
                {
                    *(*(gain+j)+i) = 0;
                }
            }
            KalmanGain(Pminus, Cr, v, K, P, gain);                         /* Kalman Gain */
            MatrixMult(Cr, *P, *K, xminus, 1, matp1);
            for(i=0; i<*P; i++)
            {
                *(*(matp1+i)) = *(*(yr+i)) - *(*(matp1+i));
            }
            MatrixMult(gain, *K, *P, matp1, 1, matk1);
            for(j=0; j<*K; j++)
            {
                *(*(filter+j)) = *(*(xminus+j)) + *(*(matk1+j));        /* Filter (x.k) */
            }
            MatrixMult(gain, *K, *P, Cr, *K, matkk);
            for(j=0; j<*K; j++)
            {
                for(jj=0; jj<*K; jj++)
                {
                    if(j != jj)
                    {
                        *(*(matkk2+j)+jj) = 0-(*(*(matkk+j)+jj));
                    }
                    if(j == jj)
                    {
                        *(*(matkk2+j)+jj) = 1-(*(*(matkk+j)+jj));
                    }
                }
            }
            MatrixMult(matkk2, *K, *K, Pminus, *K, Pk);                  /* Pk */
        }

        if(t > 0)                                                       /* All other time points */
        {
            for(j=0; j<*K; j++)
            {
                *(*(xtemp+j)) = *(*(filter+j)+(t-1));
            }
            for(i=0; i<*P; i++)
            {
                *(*(ytemp+i)) = *(*(yr+i)+(t-1));
            }
            MatrixMult(Ar, *K, *K, xtemp, 1, matk1);
            MatrixMult(Br, *K, *P, ytemp, 1, matk1b);
            for(j=0; j<*K; j++)
            {
                *(*(xminus+j)+t) = *(*(matk1+j)) + *(*(matk1b+j));      /* x.minus */
            }
            MatrixMult(Ar, *K, *K, Pk, *K, matkk);
            MatrixTrans(Ar, Art, K, K);
            MatrixMult(matkk, *K, *K, Art, *K, Pminus);
            for(j=0; j<*K; j++)
            {
                *(*(Pminus+j)+j) += 1;                                  /* P.minus */
            }

            /* Set gain = 0 */
            for(j=0; j<*K; j++)
            {
                for(i=0; i<*P; i++)
                {
                    *(*(gain+j)+i) = 0;
                }
            }
            KalmanGain(Pminus, Cr, v, K, P, gain);                         /* Kalman gain */

            for(j=0; j<*K; j++)
            {
                *(*(xtemp2+j)) = *(*(xminus+j)+t);
            }
            for(i=0; i<*P; i++)
            {
                *(*(ytemp2+i)) = *(*(yr+i)+t);
            }
            MatrixMult(Cr, *P, *K, xtemp2, 1, matp1);
            MatrixMult(Dr, *P, *P, ytemp, 1, matp1b);
            for(i=0; i<*P; i++)
            {
                *(*(ytemp+i)) = *(*(ytemp2+i)) - *(*(matp1+i)) - *(*(matp1b+i));
            }
            MatrixMult(gain, *K, *P, ytemp, 1, matk1);
            for(j=0; j<*K; j++)
            {
                *(*(filter+j)+t) = *(*(xminus+j)+t) + *(*(matk1+j));    /* Filter (x.k) */
            }
            MatrixMult(gain, *K, *P, Cr, *K, matkk);
            for(j=0; j<*K; j++)
            {
                for(jj=0; jj<*K; jj++)
                {
                    if(j != jj)
                    {
                        *(*(matkk2+j)+jj) = 0-(*(*(matkk+j)+jj));
                    }
                    if(j == jj)
                    {
                        *(*(matkk2+j)+jj) = 1-(*(*(matkk+j)+jj));
                    }
                }
            }
            MatrixMult(matkk2, *K, *K, Pminus, *K, Pk);                  /* Pk */
        }
    }

    /* Release memory */

    for(j=0; j<*K; j++)
    {
        free(*(gain+j));
        free(*(matk1+j));
        free(*(matk1b+j));
        free(*(matkk+j));
        free(*(xtemp+j));
        free(*(xtemp2+j));
        free(*(Art+j));
        free(*(matkk2+j));
    }
    for(i=0; i<*P; i++)
    {
        free(*(matp1+i));
        free(*(matp1b+i));
        free(*(ytemp+i));
        free(*(ytemp2+i));
    }
    free(Art);
    free(xtemp);
    free(xtemp2);
    free(ytemp);
    free(ytemp2);
    free(matkk);
    free(matk1);
    free(matk1b);
    free(matp1);
    free(matp1b);
    free(gain);
    free(matkk2);

}


/* Kalman Gain calculator */
void KalmanGain(double **Pminus, double **Cr, double *v, int *K, int *P, double **gain)
{
    double **matpp, **matpk, **matkp, **Crt, **inv, *det;
    int j, i;

    /* Allocate memory */
    matpp = (double**) calloc(*P, sizeof(double*));
    matpk = (double**) calloc(*P, sizeof(double*));
    matkp = (double**) calloc(*K, sizeof(double*));
    Crt = (double**) calloc(*K, sizeof(double*));
    inv = (double**) calloc(*P, sizeof(double*));
    det = (double*) calloc(1, sizeof(double));
    for(j=0; j<*K; j++)
    {
        *(matkp+j) = (double*) calloc(*P, sizeof(double));
        *(Crt+j) = (double*) calloc(*P, sizeof(double));
    }
    for(i=0; i<*P; i++)
    {
        *(matpp+i) = (double*) calloc(*P, sizeof(double));
        *(matpk+i) = (double*) calloc(*P, sizeof(double));
        *(inv+i) = (double*) calloc(*P, sizeof(double));
    }

    MatrixMult(Cr, *P, *K, Pminus, *K, matpk);
    MatrixTrans(Cr, Crt, P, K);
    MatrixMult(matpk, *P, *K, Crt, *P, matpp);

    for(i=0; i<*P; i++)
    {
        *(*(matpp+i)+i) += 1/(*(v+i));
    }
    MatrixInv(matpp, *P, inv, det);
    MatrixMult(Pminus, *K, *K, Crt, *P, matkp);
    MatrixMult(matkp, *K, *P, inv, *P, gain);

    /* Free memory */
    for(i=0; i<*P; i++)
    {
        free(*(matpp+i));
        free(*(matpk+i));
        free(*(inv+i));
    }
    for(j=0; j<*K; j++)
    {
        free(*(matkp+j));
        free(*(Crt+j));
    }
    free(matkp);
    free(matpp);
    free(matpk);
    free(Crt);
    free(inv);
    free(det);
}


/* Kalman Smoother function */
void KalmanSmoother(double **Ar, double **xminus, double **filter,
    double **Pminus, double **Pk, int *K, int *T, double **smoother,
    double **Ps)
{
    int t, j, jj;
    double **smooth, **xtemp, **matk1, **Art, **matkk;

    /* Load subprograms */
    void KalmanSmooth(double**, double**, double**, int*, double**);

    /* Allocate memory */
    smooth = (double**) calloc(*K, sizeof(double*));
    xtemp = (double**) calloc(*K, sizeof(double*));
    matk1 = (double**) calloc(*K, sizeof(double*));
    Art = (double**) calloc(*K, sizeof(double*));
    matkk = (double**) calloc(*K, sizeof(double*));
    for(j=0; j<*K; j++)
    {
        *(smooth+j) = (double*) calloc(*K, sizeof(double));
        *(xtemp+j) = (double*) calloc(1, sizeof(double));
        *(matk1+j) = (double*) calloc(1, sizeof(double));
        *(Art+j) = (double*) calloc(*K, sizeof(double));
        *(matkk+j) = (double*) calloc(*K, sizeof(double));
    }

    for(t=((*T)-1); t>=0; t--)
    {
        if(t==((*T)-1))                                                    /* Last time point */
        {
            for(j=0; j<*K; j++)
            {
                *(*(smoother+j)+t) = *(*(filter+j)+t);                     /* Smoother */
                for(jj=0; jj<*K; jj++)
                {
                    *(*(Ps+j)+jj) = *(*(Pk+j)+jj);
                }
            }
        }

        if(t<((*T)-1))                                                      /* Remaining time points */
        {
            KalmanSmooth(Pminus, Pk, Ar, K, smooth);                        /* Smooth (A.s) */
            for(j=0; j<*K; j++)
            {
                *(*(xtemp+j)) = *(*(smoother+j)+(t+1)) - *(*(xminus+j)+(t+1));
            }
            MatrixMult(smooth, *K, *K, xtemp, 1, matk1);
            for(j=0; j<*K; j++)
            {
                *(*(smoother+j)+t) = *(*(filter+j)+t) + *(*(matk1+j));      /* Smoother */
            }
            MatrixTrans(Ar, Art, K, K);
            for(j=0; j<*K; j++)
            {
                for(jj=0; jj<*K; jj++)
                {
                    *(*(matkk+j)+jj) = *(*(Ps+j)+jj) - *(*(Pminus+j)+jj);
                }
            }
            MatrixMult(smooth, *K, *K, matkk, *K, matkk);
            MatrixMult(matkk, *K, *K, Art, *K, matkk);
            for(j=0; j<*K; j++)
            {
                for(jj=0; jj<*K; jj++)
                {
                    *(*(Ps+j)+jj) = *(*(Pk+j)+jj) + *(*(matkk+j)+jj);          /* Pa */
                }
            }
        }
    }

   /* Release memory */
    for(j=0; j<*K; j++)
    {
       free(*(smooth+j));
       free(*(xtemp+j));
       free(*(matk1+j));
       free(*(Art+j));
       free(*(matkk+j));
    }
    free(smooth);
    free(xtemp);
    free(matk1);
    free(Art);
    free(matkk);

}


/* Kalman Smooth calculator */
void KalmanSmooth(double **Pminus, double **Pk, double **Ar, int *K, double **smooth)
{
    int j;
    double **Art, **matkk, **inv, *det;

    /* Allocate memory */
    Art = (double**) calloc(*K, sizeof(double*));
    matkk = (double**) calloc(*K, sizeof(double*));
    inv = (double**) calloc(*K, sizeof(double*));
    det = (double*) calloc(1, sizeof(double));
    for(j=0; j<*K; j++)
    {
        *(Art+j) = (double*) calloc(*K, sizeof(double));
        *(matkk+j) = (double*) calloc(*K, sizeof(double));
        *(inv+j) = (double*) calloc(*K, sizeof(double));
    }

    MatrixTrans(Ar, Art, K, K);
    MatrixInv(Pminus, *K, inv, det);
    MatrixMult(Pk, *K, *K, Art, *K, matkk);
    MatrixMult(matkk, *K, *K, inv, *K, smooth);

    /* Release memory */
    for(j=0; j<*K; j++)
    {
        free(*(matkk+j));
        free(*(Art+j));
        free(*(inv+j));
    }
    free(Art);
    free(matkk);
    free(inv);
    free(det);
}



void EmTypeConv(double *alpha, double *beta, double *gamma,
    double *delta, double *v, double ***x, double ***y,
    int *K, int *P, int *T, int *R, double *conv1, double *conv2,
    int *subiterations)
{
    double ***A, ***B, ***C, ***D, ***Dvar;
    int r, i, j, *KK;

    /* Load subprograms */
    void HyperMax(double*, double*, double*, double*,
        double*, double***, double***, int*, int*, int*,
        int*, double*, int*);
    void PostMeanR(double*, double*, double*, double*, double*,
        double***, double***, int*, int*, int*,	int*, double***,
        double***, double***, double***, double***);
    void VarMaxR(double***, double***, double***, double***, int*,
        int*, int*, int*, double*);

    /* Allocate space for vnew, A, B, C, D, and Dvar */
    KK = (int*) calloc(1, sizeof(int));
    *KK = 1;
    if(*K > 0) {
        *KK = *K;
    }
    A = (double***) calloc(*R, sizeof(double**));
    B = (double***) calloc(*R, sizeof(double**));
    C = (double***) calloc(*R, sizeof(double**));
    D = (double***) calloc(*R, sizeof(double**));
    Dvar = (double***) calloc(*R, sizeof(double**));
    for(r=0; r<*R; r++)
    {
        *(A+r) = (double**) calloc(*KK, sizeof(double*));
        *(B+r) = (double**) calloc(*KK, sizeof(double*));
        *(C+r) = (double**) calloc(*P, sizeof(double*));
        *(D+r) = (double**) calloc(*P, sizeof(double*));
        *(Dvar+r) = (double**) calloc(*P, sizeof(double*));
        for(j=0; j<*KK; j++)
        {
            *(*(A+r)+j) = (double*) calloc(*KK, sizeof(double));
            *(*(B+r)+j) = (double*) calloc(*P, sizeof(double));
        }
        for(i=0; i<*P; i++)
        {
            *(*(C+r)+i) = (double*) calloc(*KK, sizeof(double));
            *(*(D+r)+i) = (double*) calloc(*P, sizeof(double));
            *(*(Dvar+r)+i) = (double*) calloc(*P, sizeof(double));
        }
    }

    /* Initial estimation of hyperparameter estimates
       New estimates are put back into alpha, beta, gamma, delta
       Use convergence criterion conv1 */
    HyperMax(alpha, beta, gamma, delta, v, x, y, K, P, T, R, conv1, subiterations);
    /*Rprintf("\n");*/
    /* Estimation of gene precisions, v
       New estimate of v is put back into v */
    PostMeanR(alpha, beta, gamma, delta, v, x, y, K, P, T, R, A, B, C,
        D, Dvar);
    VarMaxR(x, y, C, D, P, R, T, K, v);
    /* Rprintf("V estimation \n"); */

    /* Final estimation of hyperparameter estimates
       New esimates are put back into alpha, beta, gamma, delta
       Use convergence criterion conv2 */
    HyperMax(alpha, beta, gamma, delta, v, x, y, K, P, T, R, conv2, subiterations);
    /*Rprintf("\n");*/
 /*  Rprintf("%f %f %f \n", *(alpha), *(alpha+1), *(alpha+2)); */

    /* Release memory */
    for(r=0; r<*R; r++)
    {
        for(j=0; j<*K; j++)
        {
            free(*(*(A+r)+j));
            free(*(*(B+r)+j));
        }
        for(i=0; i<*P; i++)
        {
            free(*(*(C+r)+i));
            free(*(*(D+r)+i));
            free(*(*(Dvar+r)+i));
        }
        free(*(A+r));
        free(*(B+r));
        free(*(C+r));
        free(*(D+r));
        free(*(Dvar+r));
    }
    free(A);
    free(B);
    free(C);
    free(D);
    free(Dvar);
    free(KK);
}



/******************************************************************/
/*  Function to conduct EM algorithm to estimate alpha - delta.   */
/*  Corresponds to the hyper.max function in R.  Inputs are       */
/*  initial values of alpha - delta, current value of v and x,    */
/*  y, K, P, T, R, and a convergence criterion.  Updated values   */
/*  of alpha - delta are returned in the original pointers.       */
/******************************************************************/


void HyperMax(double *alpha, double *beta, double *gamma, double *delta,
    double *v, double ***x, double ***y, int *K, int *P, int *T,
    int *R, double *conv, int *subiterations)
{
    int r, j, i, ii, jj, t, iter, *all, *Rchoice, *KK;
    double convergence, K2, P2;
    double **MNLNinv, **LNMNinv, **JGFGinv, **FGJGinv,
        ***HNLS, ***SNMH, ***EGFQ, ***QGJE, **matk1, **matp1,
        *tempa, *tempb, *tempc, *tempd,
        **alphanewmat, **betanewmat, **gammanewmat, **deltanewmat,
        *alphanew, *betanew, *gammanew, *deltanew,
        alphadiff, betadiff, gammadiff, deltadiff, sumnum, sumden,
        ***A, ***B, ***C, ***D, ***Dvar, *temp, *temp2;

    /* Load subprograms */
    double VecMedian(double*, int*);
    void PostMeanR(double*, double*, double*, double*, double*, double***,
        double***, int*, int*, int*, int*, double***, double***, double***,
        double***, double***);
    void SimplifyX(double*, double*, double*, double*, double*, double***,
        double***, int*, int*, int*, int*, int*,double**, double**, double**,
        double**, double***, double***, double***, double ***);

    all = (int*) calloc(1, sizeof(int));
    Rchoice = (int*) calloc(1, sizeof(int));
    KK = (int*) calloc(1, sizeof(int));
    *KK = 1;
    if(*K > 0) {
        *KK = *K;
    }
    *all = 0;
    iter = 1;
    convergence = *conv + 1.0;
    K2 = (double) (*K)/2.0;
    P2 = (double) (*P)/2.0;

    /* Allocate MNLNinv, LNMNinv, JGFGinv, FGJGinv, HNLS, SNMH, EGFQ, and QGJE and set to 0 */
    MNLNinv = (double**) calloc(*KK, sizeof(double*));
    LNMNinv = (double**) calloc(*P, sizeof(double*));
    JGFGinv = (double**) calloc(*KK, sizeof(double*));
    FGJGinv = (double**) calloc(*P, sizeof(double*));
    HNLS = (double***) calloc(*KK, sizeof(double**));
    SNMH = (double***) calloc(*KK, sizeof(double**));
    EGFQ = (double***) calloc(*P, sizeof(double**));
    QGJE = (double***) calloc(*P, sizeof(double**));

    matk1 = (double**) calloc(*KK, sizeof(double*));
    matp1 = (double**) calloc(*P, sizeof(double*));
    tempa = (double*) calloc(*KK, sizeof(double));
    tempb = (double*) calloc(*P, sizeof(double));
    tempc = (double*) calloc(*KK, sizeof(double));
    tempd = (double*) calloc(*P, sizeof(double));

    alphanewmat = (double**) calloc(*KK, sizeof(double*));
    betanewmat = (double**) calloc(*P, sizeof(double*));
    gammanewmat = (double**) calloc(*KK, sizeof(double*));
    deltanewmat = (double**) calloc(*P, sizeof(double*));
    alphanew = (double*) calloc(*KK, sizeof(double));
    betanew = (double*) calloc(*P, sizeof(double));
    gammanew = (double*) calloc(*KK, sizeof(double));
    deltanew = (double*) calloc(*P, sizeof(double));

    for(j=0; j<*KK; j++)
    {
        *(MNLNinv+j) = (double*) calloc(*KK, sizeof(double));
        *(JGFGinv+j) = (double*) calloc(*KK, sizeof(double));
        *(HNLS+j) = (double**) calloc(*KK, sizeof(double*));
        *(SNMH+j) = (double**) calloc(*P, sizeof(double*));
        for(jj=0; jj<*KK; jj++)
        {
            *(*(HNLS+j)+jj) = (double*) calloc(1, sizeof(double));
        }
        for(i=0; i<*P; i++)
        {
            *(*(SNMH+j)+i) = (double*) calloc(1, sizeof(double));
        }
        *(matk1+j) = (double*) calloc(1, sizeof(double));
        *(alphanewmat+j) = (double*) calloc(*R, sizeof(double));
        *(gammanewmat+j) = (double*) calloc(*R, sizeof(double));
    }
    for(i=0; i<*P; i++)
    {
        *(LNMNinv+i) = (double*) calloc(*P, sizeof(double));
        *(FGJGinv+i) = (double*) calloc(*P, sizeof(double));
        *(EGFQ+i) = (double**) calloc(*KK, sizeof(double*));
        for(j=0; j<*KK; j++)
        {
            *(*(EGFQ+i)+j) = (double*) calloc(1, sizeof(double));
        }
        *(QGJE+i) = (double**) calloc(*P, sizeof(double*));
        for(ii=0; ii<*P; ii++)
        {
            *(*(QGJE+i)+ii) = (double*) calloc(1, sizeof(double));
        }
        *(matp1+i) = (double*) calloc(1, sizeof(double));
        *(betanewmat+i) = (double*) calloc(*R, sizeof(double));
        *(deltanewmat+i) = (double*) calloc(*R, sizeof(double));
    }

    /* Rprintf("Hypermax loop! "); */

    /* If estimating x's */
    if(*K > 0)
    {
        while(convergence > *conv)
        {
            if(iter > *subiterations) break;

            for(r=0; r<*R; r++)
            {
                *Rchoice = r;
                /* Set MNLNinv, LNMNinv, JGFGinv, FGJGinv, HNLS, SNMH, EGFQ, QGJE to 0 */
                for(j=0; j<*K; j++)
                {
                    for(jj=0; jj<*K; jj++)
                    {
                        *(*(MNLNinv+j)+jj)=0;
                        *(*(JGFGinv+j)+jj)=0;
                        *(*(*(HNLS+j)+jj))=0;
                    }
                    for(i=0; i<*P; i++)
                    {
                        *(*(*(SNMH+j)+i))=0;
                    }
                    *(*(matk1+j)) = 0;
                }
                for(i=0; i<*P; i++)
                {
                    for(ii=0; ii<*P; ii++)
                    {
                        *(*(LNMNinv+i)+ii)=0;
                        *(*(FGJGinv+i)+ii)=0;
                        *(*(*(QGJE+i)+ii))=0;
                    }
                    for(j=0; j<*K; j++)
                    {
                        *(*(*(EGFQ+i)+j))=0;
                    }
                }

                SimplifyX(alpha, beta, gamma, delta, v, x, y, K, P,             /* Simplify function */
                    T, all, Rchoice, MNLNinv, LNMNinv, JGFGinv,
                    FGJGinv, HNLS, SNMH, EGFQ, QGJE);

                for(j=0; j<*K; j++)
                {
                    *(tempa+j) = 0;
                    *(tempc+j) = 0;
                }
                for(i=0; i<*P; i++)
                {
                    *(tempb+i) = 0;
                    *(tempd+i) = 0;
                }

                for(j=0; j<*K; j++)
                {
                      MatrixMult(MNLNinv, *K, *K, *(HNLS+j), 1, matk1);
                    for(jj=0; jj<*K; jj++)
                    {
                        *(tempa+jj) += (*(*(matk1+jj))) * (*(*(matk1+jj)));     /* temp.a */
                    }
                    MatrixMult(LNMNinv, *P, *P, *(SNMH+j), 1, matp1);
                    for(i=0; i<*P; i++)
                    {
                        *(tempb+i) += (*(*(matp1+i))) * (*(*(matp1+i)));        /* temp.b */
                    }
                }

                for(i=0; i<*P; i++)
                {
                    MatrixMult(JGFGinv, *K, *K, *(EGFQ+i), 1, matk1);
                    for(j=0; j<*K; j++)
                    {
                        *(*(matk1+j)) = (*(*(matk1+j))) * (*(*(matk1+j))) * (*(v+i));
                        *(tempc+j) += (*(*(matk1+j)));                          /* temp.c */
                    }
                    MatrixMult(FGJGinv, *P, *P, *(QGJE+i), 1, matp1);
                    for(ii=0; ii<*P; ii++)
                    {
                        *(*(matp1+ii)) = (*(*(matp1+ii))) * (*(*(matp1+ii))) * (*(v+i));
                        *(tempd+ii) += (*(*(matp1+ii)));                        /* temp.d */
                    }
                }

                for(j=0; j<*K; j++)
                {
                    *(*(alphanewmat+j)+r) = K2 * (1.0/(K2 * (*(*(MNLNinv+j)+j)) + (1.0/2.0)*(*(tempa+j))));
                    *(*(gammanewmat+j)+r) = P2 * (1.0/(P2 * (*(*(JGFGinv+j)+j)) + (1.0/2.0)*(*(tempc+j))));
                }
                for(i=0; i<*P; i++)
                {
                    *(*(betanewmat+i)+r) = K2 * (1/(K2 * (*(*(LNMNinv+i)+i)) + (1.0/2.0)*(*(tempb+i))));
                    *(*(deltanewmat+i)+r) = P2 * (1/(P2 * (*(*(FGJGinv+i)+i)) + (1.0/2.0)*(*(tempd+i))));
                }

                /* Rprintf("R%d ", r); */
            }
            /*Rprintf("\n");*/

            for(j=0; j<*K; j++)
            {
                *(alphanew+j) = VecMedian(*(alphanewmat+j), R);             /* alpha.new */
                *(gammanew+j) = VecMedian(*(gammanewmat+j), R);             /* beta.new */
            }
            for(i=0; i<*P; i++)
            {
                *(betanew+i) = VecMedian(*(betanewmat+i),R);                /* gamma.new */
                *(deltanew+i) = VecMedian(*(deltanewmat+i),R);              /* delta.new */
            }

            /* Check convergence */

            sumnum = 0;
            sumden = 0;
            for(j=0; j<*K; j++)
            {
                sumden += (*(alpha+j)) * (*(alpha+j));
            }
            for(j=0; j<*K; j++)
            {
                sumnum += ((*(alpha+j) - *(alphanew+j))*(*(alpha+j) - *(alphanew+j)))/sumden;
            }
            alphadiff = sqrt(sumnum);                                       /* alpha.diff */

            sumnum = 0;
            sumden = 0;
            for(i=0; i<*P; i++)
            {
                sumden += (*(beta+i)) * (*(beta+i));
            }
            for(i=0; i<*P; i++)
            {
                sumnum += ((*(beta+i) - *(betanew+i))*(*(beta+i) - *(betanew+i)))/sumden;
            }
            betadiff = sqrt(sumnum);                                       /* beta.diff */

            sumnum = 0;
            sumden = 0;
            for(j=0; j<*K; j++)
            {
                sumden += (*(gamma+j)) * (*(gamma+j));
            }
            for(j=0; j<*K; j++)
            {
                sumnum += ((*(gamma+j) - *(gammanew+j))*(*(gamma+j) - *(gammanew+j)))/sumden;
            }
            gammadiff = sqrt(sumnum);                                       /* gamma.diff */

            sumnum = 0;
            sumden = 0;
            for(i=0; i<*P; i++)
            {
                sumden += (*(delta+i)) * (*(delta+i));
            }
            for(i=0; i<*P; i++)
            {
                sumnum += ((*(delta+i) - *(deltanew+i))*(*(delta+i) - *(deltanew+i)))/sumden;
            }
            deltadiff = sqrt(sumnum);                                       /* delta.diff */

            /* Find maximum convergence criterion */
            convergence = alphadiff;                                        /* Convergence criterion */
            if(betadiff>convergence) {convergence = betadiff;}
            if(gammadiff>convergence) {convergence = gammadiff;}
            if(deltadiff>convergence) {convergence = deltadiff;}

            /* Rprintf("Check: "); */

            for(j=0; j<*K; j++)                                             /* Update hyperparameters */
            {
                *(alpha+j) = *(alphanew+j);
                *(gamma+j) = *(gammanew+j);
            }
            for(i=0; i<*P; i++)
            {
                *(beta+i) = *(betanew+i);
                *(delta+i) = *(deltanew+i);
            }
            /*if(*K > 0)
            {
                Rprintf("%d ", iter);
            }*/
            iter++;                                                            /* Update iteration number */
        }
        /*Rprintf("\n");*/
    }
    /* Free MNLNinv, LNMNinv, JGFGinv, FGJGinv, HNLS, SNMH, EGFQ, and QGJE */

    for(j=0; j<*KK; j++)
    {
        for(jj=0; jj<*KK; jj++)
        {
            free(*(*(HNLS+j)+jj));
        }
        for(i=0; i<*P; i++)
        {
            free(*(*(SNMH+j)+i));
        }
        free(*(HNLS+j));
        free(*(SNMH+j));
        free(*(MNLNinv+j));
        free(*(JGFGinv+j));
        free(*(matk1+j));
        free(*(alphanewmat+j));
        free(*(gammanewmat+j));
    }

    for(i=0; i<*P; i++)
    {
        for(j=0; j<*KK; j++)
        {
            free(*(*(EGFQ+i)+j));
        }
        for(ii=0; ii<*P; ii++)
        {
            free(*(*(QGJE+i)+ii));
        }
        free(*(EGFQ+i));
        free(*(QGJE+i));
        free(*(LNMNinv+i));
        free(*(FGJGinv+i));
        free(*(matp1+i));
        free(*(betanewmat+i));
        free(*(deltanewmat+i));
    }
    free(HNLS);
    free(SNMH);
    free(MNLNinv);
    free(JGFGinv);
    free(EGFQ);
    free(QGJE);
    free(LNMNinv);
    free(FGJGinv);

    free(matk1);
    free(matp1);
    free(tempa);
    free(tempb);
    free(tempc);
    free(tempd);
    free(alphanewmat);
    free(betanewmat);
    free(gammanewmat);
    free(deltanewmat);
    free(alphanew);
    free(betanew);
    free(gammanew);
    free(deltanew);


    /* If NOT estimating x's */
    A = (double***) calloc(*R, sizeof(double**));
    B = (double***) calloc(*R, sizeof(double**));
    C = (double***) calloc(*R, sizeof(double**));
    D = (double***) calloc(*R, sizeof(double**));
    Dvar = (double***) calloc(*R, sizeof(double**));
    temp = (double*) calloc(*P, sizeof(double));
    temp2 = (double*) calloc(*P, sizeof(double));
    deltanewmat = (double**) calloc(*P, sizeof(double*));
    deltanew = (double*) calloc(*P, sizeof(double));
    matp1 = (double**) calloc(*P, sizeof(double*));

    for(i=0; i<*P; i++)
    {
        *(deltanewmat+i) = (double*) calloc(*R, sizeof(double));
        *(matp1+i) = (double*) calloc(1, sizeof(double));
    }
    for(r=0; r<*R; r++)
    {
        *(A+r) = (double**) calloc(*KK, sizeof(double*));
        *(B+r) = (double**) calloc(*KK, sizeof(double*));
        *(C+r) = (double**) calloc(*P, sizeof(double*));
        *(D+r) = (double**) calloc(*P, sizeof(double*));
        *(Dvar+r) = (double**) calloc(*P, sizeof(double*));
        for(j=0; j<*KK; j++)
        {
            *(*(A+r)+j) = (double*) calloc(*KK, sizeof(double));
            *(*(B+r)+j) = (double*) calloc(*P, sizeof(double));
        }
        for(i=0; i<*P; i++)
        {
            *(*(C+r)+i) = (double*) calloc(*KK, sizeof(double));
            *(*(D+r)+i) = (double*) calloc(*P, sizeof(double));
            *(*(Dvar+r)+i) = (double*) calloc(*P, sizeof(double));
        }
    }

    if(*K == 0)
    {
        while(convergence > *conv)
        {
            if(iter > *subiterations) break;

            PostMeanR(alpha, beta, gamma, delta, v, x, y, K, P, T, R,
                A, B, C, D, Dvar);
            for(r=0; r<*R; r++)
            {
                for(i=0; i<*P; i++)
                {
                    *(temp+i) = 0;
                }
                for(i=0; i<*P; i++)
                {
                    for(ii=0; ii<*P; ii++)
                    {
                        *(temp2+ii) = 0;
                    }
                    for(t=0; t<((*T)-1); t++)
                    {
                        for(ii=0; ii<*P; ii++)
                        {
                            *(temp2+ii) += (*(*(*(y+r)+ii)+t)) * (*(*(*(y+r)+i)+(t+1)));
                        }
                    }
                    /* Set matp1 = temp2 */
                    for(ii=0; ii<*P; ii++)
                    {
                        *(*(matp1+ii)) = *(temp2+ii);
                    }

                    MatrixMult(*(Dvar+r), *P, *P, matp1, 1, matp1);
                    for(ii=0; ii<*P; ii++)
                    {
                        *(temp+ii) += (*(v+i)) * (*(*(matp1+ii))) * (*(*(matp1+ii)));
                    }
                }

                for(i=0; i<*P; i++)
                {
                    *(*(deltanewmat+i)+r) = P2 * (1/(P2* (*(*(*(Dvar+r)+i)+i)) + (1.0/2.0) * (*(temp+i))));
                }
            }
            for(i=0; i<*P; i++)
            {
                *(deltanew+i) = VecMedian(*(deltanewmat+i),R);                  /* delta.new */
            }

            sumnum = 0;
            sumden = 0;
            for(i=0; i<*P; i++)
            {
                sumden += (*(delta+i)) * (*(delta+i));
            }
            for(i=0; i<*P; i++)
            {
                sumnum += ((*(delta+i) - *(deltanew+i))*(*(delta+i) - *(deltanew+i)))/sumden;
            }
            convergence = sqrt(sumnum);                                         /* Convergence criterion */
 /*           Rprintf("%f, ", convergence); */
            for(i=0; i<*P; i++)                                                 /* Update delta */
            {
                *(delta+i) = *(deltanew+i);
            }
    /*        Rprintf("C iteration %d \n", iter); */
            iter++;                                                             /* Update iteration number */
        }
    }
    /* Release memory */
    for(r=0; r<*R; r++)
    {
        for(j=0; j<*KK; j++)
        {
            free(*(*(A+r)+j));
            free(*(*(B+r)+j));
        }
        for(i=0; i<*P; i++)
        {
            free(*(*(C+r)+i));
            free(*(*(D+r)+i));
            free(*(*(Dvar+r)+i));
        }
        free(*(A+r));
        free(*(B+r));
        free(*(C+r));
        free(*(D+r));
        free(*(Dvar+r));
    }
    for(i=0; i<*P; i++)
    {
        free(*(deltanewmat+i));
        free(*(matp1+i));
    }
    free(A);
    free(B);
    free(C);
    free(D);
    free(Dvar);
    free(temp);
    free(temp2);
    free(deltanewmat);
    free(deltanew);
    free(matp1);
    free(KK);
}


/******************************************************************/
/*  This is the counterpart of Posterior_Mean.R.  It takes as     */
/*  inputs alpha, beta, gamma, delta, v, x, y, K, P, T, and R     */
/*  and returns the value of the posterior mean of A, B, C, and D */
/*  (in the case where x's are estimated), or the posterior mean  */
/*  and variance of D (in the case where x's are not estimated).  */
/******************************************************************/


void PostMeanR(double *alpha, double *beta, double *gamma, double *delta,
	double *v, double ***x, double ***y, int *K, int *P, int *T,
	int *R,	double ***A, double ***B, double ***C, double ***D, double ***Dvar)
{
	int r, i, ii, *all, j, jj, *Rchoice, *KK;
	double **Dmean, **Dv;
    double **MNLNinv, **LNMNinv, **JGFGinv, **FGJGinv, ***HNLS,
            ***SNMH, ***EGFQ, ***QGJE, **matk1, **matp1;

	/* Load subprograms */
	void SimplifyNoX(double*, double*, double***, int*, int*, int*, int*,
        	double**, double**);
	void SimplifyX(double*, double*, double*, double*, double*, double***,
		double***, int*, int*, int*, int*, int*, double**, double**, double**,
        	double**, double***, double***, double***,double***);

	all = (int*) calloc(1, sizeof(int));
	Rchoice = (int*) calloc(1, sizeof(int));
	KK = (int*) calloc(1, sizeof(int));
	*all = 0;
	*KK = 1;
    if(*K > 0)
    {
        *KK = *K;
    }

    /* Allocate memory */
    Dmean = (double**) calloc(*P, sizeof(double*));
    Dv = (double**) calloc(*P, sizeof(double*));
    for(i=0; i<*P; i++)
    {
        *(Dmean + i) = (double*) calloc(*P, sizeof(double));
        *(Dv + i) = (double*) calloc(*P, sizeof(double));
    }
    /* Allocate MNLNinf, LNMNinf, JGFGinv, FGJGinv, HNLS, SNMH, EGFQ, QFJE */
    MNLNinv = (double**) calloc(*KK, sizeof(double*));
    LNMNinv = (double**) calloc(*P, sizeof(double*));
    JGFGinv = (double**) calloc(*KK, sizeof(double*));
    FGJGinv = (double**) calloc(*P, sizeof(double*));
    HNLS = (double***) calloc(*KK, sizeof(double**));
    SNMH = (double***) calloc(*KK, sizeof(double**));
    EGFQ = (double***) calloc(*P, sizeof(double**));
    for(j=0; j<*KK; j++)
    {
        *(MNLNinv+j) = (double*) calloc(*KK, sizeof(double));
        *(JGFGinv+j) = (double*) calloc(*KK, sizeof(double));
        *(HNLS+j) = (double**) calloc(*KK, sizeof(double*));
        *(SNMH+j) = (double**) calloc(*P, sizeof(double*));
        for(jj=0; jj<*KK; jj++)
        {
            *(*(HNLS+j)+jj) = (double*) calloc(1, sizeof(double));
        }
        for(i=0; i<*P; i++)
        {
            *(*(SNMH+j)+i) = (double*) calloc(1, sizeof(double));
        }
    }
    for(i=0; i<*P; i++)
    {
        *(LNMNinv+i) = (double*) calloc(*P, sizeof(double));
        *(FGJGinv+i) = (double*) calloc(*P, sizeof(double));
        *(EGFQ+i) = (double**) calloc(*KK, sizeof(double*));
        for(j=0; j<*KK; j++)
        {
            *(*(EGFQ+i)+j) = (double*) calloc(1, sizeof(double));
        }
    }
    QGJE = (double***) calloc(*P, sizeof(double**));
    for(i=0; i<*P; i++)
    {
        *(QGJE+i) = (double**) calloc(*P, sizeof(double*));
        for(ii=0; ii<*P; ii++)
        {
            *(*(QGJE+i)+ii) = (double*) calloc(1, sizeof(double));
        }
    }

    matk1 = (double**) calloc(*KK, sizeof(double*));
    matp1 = (double**) calloc(*P, sizeof(double*));
    for(j=0; j<*KK; j++)
    {
        *(matk1+j) = (double*) calloc(1, sizeof(double));
    }
    for(i=0; i<*P; i++)
    {
        *(matp1+i) = (double*) calloc(1, sizeof(double));
    }


	/* Begin function -- no x's */


	if(*K==0)
	{
		for(r=0; r<*R; r++)
		{
			*Rchoice = r;
			/* Set Dmean and Dv to 0 */
			for(i=0; i<*P; i++)
			{
			    for(ii=0; ii<*P; ii++)
			    {
			        *(*(Dmean+i)+ii) = 0;
			        *(*(Dv+i)+ii) = 0;
			    }
			}
			SimplifyNoX(delta, v, y, P, T, Rchoice, all, Dmean, Dv);
			for(i=0; i<*P; i++)
			{
				for(ii=0; ii<*P; ii++)
				{
					*(*(*(D+r)+i)+ii) = *(*(Dmean+i)+ii);
					*(*(*(Dvar+r)+i)+ii) = *(*(Dv+i)+ii);
				}
			}
		}
	}
	/* Begin function -- x's */
	if(*K > 0)
	{
		for(r=0; r<*R; r++)
		{
			*Rchoice = r;

            /* Set MNLNinv, LNMNinv, JGFGinv, FGJGinv, HNLS, SNMH, EGFQ, QGJE to 0 */
            for(j=0; j<*K; j++)
            {
                for(jj=0; jj<*K; jj++)
                {
                    *(*(MNLNinv+j)+jj)=0;
                    *(*(JGFGinv+j)+jj)=0;
                    *(*(*(HNLS+j)+jj))=0;
                }
                for(i=0; i<*P; i++)
                {
                    *(*(*(SNMH+j)+i))=0;
                }
            }
            for(i=0; i<*P; i++)
            {
                for(ii=0; ii<*P; ii++)
                {
                    *(*(LNMNinv+i)+ii)=0;
                    *(*(FGJGinv+i)+ii)=0;
                    *(*(*(QGJE+i)+ii))=0;
                }
                for(j=0; j<*K; j++)
                {
                    *(*(*(EGFQ+i)+j))=0;
                }
            }

			SimplifyX(alpha, beta, gamma, delta, v, x, y, K, P,
				T, all, Rchoice, MNLNinv, LNMNinv, JGFGinv,
                FGJGinv, HNLS, SNMH, EGFQ,QGJE);
			for(j=0; j<*K; j++)
			{
				MatrixMult(MNLNinv, *K, *K, *(HNLS+j), 1, matk1);
				MatrixMult(LNMNinv, *P, *P, *(SNMH+j), 1, matp1);
				for(jj=0; jj<*K; jj++)
				{
					*(*(*(A+r)+j)+jj) = *(*(matk1+jj));
				}
				for(i=0; i<*P; i++)
				{
					*(*(*(B+r)+j)+i) = *(*(matp1+i));
				}
			}
			for(i=0; i<*P; i++)
			{
				MatrixMult(JGFGinv, *K, *K, *(EGFQ+i), 1, matk1);
				MatrixMult(FGJGinv, *P, *P, *(QGJE+i), 1, matp1);
				for(j=0; j<*K; j++)
				{
					*(*(*(C+r)+i)+j) = *(*(matk1+j));
				}
				for(ii=0; ii<*P; ii++)
				{
					*(*(*(D+r)+i)+ii) = *(*(matp1+ii));
				}
			}
		}
	}

	/* Free memory */

	/* If no x's */
    for(i=0; i<*P; i++)
    {
        free(*(Dmean+i));
        free(*(Dv+i));
    }
    free(Dmean);
    free(Dv);
    for(j=0; j<*KK; j++)
    {
        for(jj=0; jj<*KK; jj++)
        {
            free(*(*(HNLS+j)+jj));
        }
        for(i=0; i<*P; i++)
        {
            free(*(*(SNMH+j)+i));
        }
        free(*(HNLS+j));
        free(*(SNMH+j));
        free(*(MNLNinv+j));
        free(*(JGFGinv+j));
    }

    for(i=0; i<*P; i++)
    {
        for(j=0; j<*KK; j++)
        {
            free(*(*(EGFQ+i)+j));
        }
        for(ii=0; ii<*P; ii++)
        {
            free(*(*(QGJE+i)+ii));
        }
        free(*(EGFQ+i));
        free(*(QGJE+i));
        free(*(LNMNinv+i));
        free(*(FGJGinv+i));
    }
    free(HNLS);
    free(SNMH);
    free(MNLNinv);
    free(JGFGinv);
    free(EGFQ);
    free(QGJE);
    free(LNMNinv);
    free(FGJGinv);

    for(j=0; j<*KK; j++)
    {
        free(*(matk1+j));
    }
    for(i=0; i<*P; i++)
    {
        free(*(matp1+i));
    }
    free(matk1);
    free(matp1);
	free(all);
	free(Rchoice);
	free(KK);
}


/* *************************************************** */
/* Function to calculate the median variance estimator */
/* This function is the C counterpart of the program   */
/* Variance_Estimator.R                                */
/* *************************************************** */

void VarMaxR(double ***x, double ***y, double ***C, double ***D, int *P,
	int *R, int *T, int *K, double *vEst)
{
	int r, t, i, j, *KK;
	double **vTemp, **vTempt, *temp, **tempDy, **tempCx, **tempyold,
		*tempyDy, *tempyCx, **tempxold, *tempyCxDy;

	/* Load subprograms */
	double VecMedian(double*, int*);

	/* Allocate temporary and vTemp matrices */
	KK = (int*) calloc(1, sizeof(int));
	*KK = 1;
	if(*K > 0)
	{
	    *KK = *K;
	}
	vTemp = (double**) calloc(*R, sizeof(double*));
	vTempt = (double**) calloc(*P, sizeof(double*));
	temp = (double*) calloc(*P, sizeof(double));
	tempDy = (double**) calloc(*P, sizeof(double*));
	tempCx = (double**) calloc(*P, sizeof(double*));
	tempyold = (double**) calloc(*P, sizeof(double*));
	tempyDy = (double*) calloc(*P, sizeof(double));
	tempyCx = (double*) calloc(*P, sizeof(double));
	tempxold = (double**) calloc(*KK, sizeof(double*));
	tempyCxDy = (double*) calloc(*P, sizeof(double));

	for(r=0; r<*R; r++)
	{
		*(vTemp+r) = (double*) calloc(*P, sizeof(double));
	}
	for(i=0; i<*P; i++)
	{
		*(tempDy+i) = (double*) calloc(1, sizeof(double));
		*(tempCx+i) = (double*) calloc(1, sizeof(double));
		*(vTempt+i) = (double*) calloc(*R, sizeof(double));
		*(tempyold+i) = (double*) calloc(1, sizeof(double));
	}
	for(j=0; j<*K; j++)
	{
		*(tempxold+j) = (double*) calloc(1, sizeof(double));
	}

	/* Calculate the variance estimator if no x's */
	if(*K == 0)
	{
		for(r=0; r<*R; r++)
		{
			for(i=0; i<*P; i++)
			{
				*(temp+i) = 0;
			}
			for(t=1; t<*T; t++)
			{
				for(i=0; i<*P; i++)
				{
					*(*(tempyold+i))=*(*(*(y+r)+i)+(t-1));
				}
				MatrixMult(*(D+r),*P,*P,tempyold,1,tempDy);
				for(i=0; i<*P; i++)
				{
					*(tempyDy+i)=*(*(*(y+r)+i)+t) - *(*(tempDy+i));
					*(temp+i) += (*(tempyDy+i)) * (*(tempyDy+i));
				}
			}

			for(i=0; i<*P; i++)
			{
				*(*(vTemp+r)+i) = (*(temp+i)) / ((*T)-1);
			}
		}
	}

	/* Calculate the variance estimator if x's */
        if(*K > 0)
        {
                for(r=0; r<*R; r++)
                {
                        for(i=0; i<*P; i++)
                        {
                                *(temp+i) = 0;
                        }
			for(j=0; j<*K; j++)
			{
				*(*(tempxold+j)) = *(*(*(x+r)+j));
			}
			MatrixMult(*(C+r),*P,*K,tempxold,1,tempCx);
			for(i=0; i<*P; i++)
			{
				*(tempyCx+i)=*(*(*(y+r)+i)) - *(*(tempCx+i));
				*(temp+i) += (*(tempyCx+i)) * (*(tempyCx+i));
			}

                        for(t=1; t<*T; t++)
                        {
                                for(i=0; i<*P; i++)
                                {
                                        *(*(tempyold+i))=*(*(*(y+r)+i)+(t-1));
                                }
				for(j=0; j<*K; j++)
				{
					*(*(tempxold+j))=*(*(*(x+r)+j)+t);
				}
                MatrixMult(*(D+r),*P,*P,tempyold,1,tempDy);
				MatrixMult(*(C+r),*P,*K,tempxold,1,tempCx);
                                for(i=0; i<*P; i++)
                                {
                                        *(tempyCxDy+i)=*(*(*(y+r)+i)+t) - *(*(tempCx+i)) - *(*(tempDy+i));
                                        *(temp+i) += (*(tempyCxDy+i)) * (*(tempyCxDy+i));
                                }
                        }

                        for(i=0; i<*P; i++)
                        {
                                *(*(vTemp+r)+i) = ( *(temp+i) / (*T) );
                        }
                }
        }

	/* Take inverse of all elements of vTemp matrix */
	for(r=0; r<*R; r++)
	{
		for(i=0; i<*P; i++)
		{
			*(*(vTemp+r)+i) = 1/(*(*(vTemp+r)+i));
		}
	}

	/* Rewrite the vTemp matrix as its transpose, find the median */
	MatrixTrans(vTemp, vTempt, R, P);
	for(i=0; i<*P; i++)
	{
		*(vEst+i) = VecMedian(*(vTempt+i), R);
	}

	/* Free memory */
	for(r=0; r<*R; r++)
	{
		free(*(vTemp+r));
	}
	for(i=0; i<*P; i++)
	{
		free(*(tempDy+i));
		free(*(tempCx+i));
		free(*(vTempt+i));
		free(*(tempyold+i));
	}
	for(j=0; j<*KK; j++)
	{
		free(*(tempxold+j));
	}
	free(tempyCx);
	free(tempxold);
	free(tempyCxDy);
	free(vTemp);
	free(vTempt);
	free(temp);
	free(tempDy);
	free(tempCx);
	free(tempyold);
	free(tempyDy);
    free(KK);
}


/*********************************************************************/
/*  Function to calculate the median of a vector of a given length.  */
/*  Will return the median as a double.                              */
/*********************************************************************/

double VecMedian(double *vector, int *length)
{
	int index1, index2;
 	double len, tempnum = 1.0, median = 0;
	len = (double) (*length)/tempnum;

	/* Sort the vector using the rsort function in R */
	R_rsort(vector, *length);

	/* Odd case */
	if((len)/2 != floor((len)/2))
	{
		index1 = (int) floor((len)/2);
		median = *(vector + index1);
	}

	/* Even case */
	if((len)/2 == floor((len)/2))
	{
		index1 = (int) (((len)/2)-1);
		index2 = (int) ((len)/2);
		median = ((*(vector + index1)) + (*(vector + index2)))/2;
	}
	return(median);
}




/*************************************************/
/*  This is the C version of Simplifying_code.R  */
/*  in the case where x's are estimated.         */
/*************************************************/

void SimplifyX(double *alpha, double *beta, double *gamma, double *delta,
	double *v, double ***x, double ***y, int *K, int *P, int *T, int *all,
	int *Rchoice,
	double **MNLNinv, double **LNMNinv, double **JGFGinv,
	double **FGJGinv, double ***HNLS, double ***SNMH, double ***EGFQ,
	double ***QGJE)
{
	int rlower = 0, rupper = 0, i, j, ii, jj, r, t, *one;
	double **xXx, **yXx, **yXy, **xXy, ***H, ***S, **M, **N, **L, **xtemp,
		**txtemp, **ytemp, **tytemp, *xsingle, *ysingle, **J, **G,
		**F, ***E, ***Q;
	double **Linv, **Minv, **Finv, **Jinv, *det, **Nt, **Gt,
		**matkk, **matpp, **matkp, **matpk;

	/* Allocate memory, set M, N, L, H, and S = 0 */
	/* Allocate memory, set J, G, F, E, and Q = 0 */

	one = (int*) calloc(1, sizeof(int));
	det = (double*) calloc(1, sizeof(double));
    H = (double***) calloc(*K, sizeof(double**));
    S = (double***) calloc(*K, sizeof(double**));
    M = (double**) calloc(*K, sizeof(double*));
    N = (double**) calloc(*P, sizeof(double*));
    L = (double**) calloc(*P, sizeof(double*));
    J = (double**) calloc(*K, sizeof(double*));
    G = (double**) calloc(*K, sizeof(double*));
    F = (double**) calloc(*P, sizeof(double*));
    Q = (double***) calloc(*P, sizeof(double**));
    E = (double***) calloc(*P, sizeof(double**));

	xXx = (double**) calloc(*K, sizeof(double*));
	yXx = (double**) calloc(*P, sizeof(double*));
	yXy = (double**) calloc(*P, sizeof(double*));
	xXy = (double**) calloc(*K, sizeof(double*));
	xtemp = (double**) calloc(*K, sizeof(double*));
	txtemp = (double**) calloc(1, sizeof(double*));
	ytemp = (double**) calloc(*P, sizeof(double*));
	tytemp = (double**) calloc(1, sizeof(double*));
	xsingle = (double*) calloc(1, sizeof(double));
	ysingle = (double*) calloc(1, sizeof(double));


	for(i=0; i<*P; i++)
	{
		*(yXx + i) = (double*) calloc(*K, sizeof(double));
		*(yXy + i) = (double*) calloc(*P, sizeof(double));
        *(ytemp + i) = (double*) calloc(1, sizeof(double));
		*(N + i) = (double*) calloc(*K, sizeof(double));
		*(L + i) = (double*) calloc(*P, sizeof(double));
		*(F + i) = (double*) calloc(*P, sizeof(double));
		*(E + i) = (double**) calloc(*K, sizeof(double*));
		*(Q + i) = (double**) calloc(*P, sizeof(double*));
		for(j=0; j<*K; j++)
		{
			*(*(N+i)+j) = 0;
			*(*(E+i)+j) = (double*) calloc(1, sizeof(double));
			*(*(*(E+i)+j)) = 0;
		}
		for(ii=0; ii<*P; ii++)
		{
			*(*(L+i)+ii) = 0;
			*(*(F+i)+ii) = 0;
			*(*(Q+i)+ii) = (double*) calloc(1, sizeof(double));
			*(*(*(Q+i)+ii)) = 0;
		}
	}
	for(j=0; j<*K; j++)
	{
		*(xXx + j) = (double*) calloc(*K, sizeof(double));
		*(xXy + j) = (double*) calloc(*P, sizeof(double));
		*(xtemp + j) = (double*) calloc(1, sizeof(double));
        *(H + j) = (double**) calloc(*K, sizeof(double*));
        *(S + j) = (double**) calloc(*P, sizeof(double*));
		*(M + j) = (double*) calloc(*K, sizeof(double));
		*(J + j) = (double*) calloc(*K, sizeof(double));
		*(G + j) = (double*) calloc(*P, sizeof(double));
		for(jj=0; jj<*K; jj++)
		{
			*(*(H+j)+jj) = (double*) calloc(1, sizeof(double));
			*(*(M+j)+jj) = 0;
			*(*(*(H+j)+jj)) = 0;
			*(*(J+j)+jj) = 0;
		}
		for(i=0; i<*P; i++)
		{
			*(*(S+j)+i) = (double*) calloc(1, sizeof(double));
			*(*(*(S+j)+i)) = 0;
			*(*(G+j)+i) = 0;
		}
	}
	*(txtemp) = (double*) calloc(*K, sizeof(double));
	*(tytemp) = (double*) calloc(*P, sizeof(double));

	/* Set limits of function */

	if(*all == 1)
	{
		rlower = 0;
		rupper = *Rchoice;
	}
	if(*all == 0)
	{
		rlower = *Rchoice;
		rupper = rlower + 1;
	}
	*one = 1;

	/* Begin function part 1 (corresponds to simplify.1 in R code) */

	for(r=rlower; r<rupper; r++)
	{
		for(t=0; t<((*T)-1); t++)
		{
			for(i=0; i<*P; i++)
			{
				*(*(ytemp+i)) = *(*(*(y+r)+i)+t);
				*(*(tytemp)+i) = *(*(*(y+r)+i)+t);
			}
			for(j=0; j<*K; j++)
            {
                *(*(xtemp+j)) = *(*(*(x+r)+j)+(t+1));
            }
            MatrixMult(xtemp,*K,1,tytemp,*P,xXy);           /* temp.kp */

  			for(j=0; j<*K; j++)
            {
                *(*(xtemp+j)) = *(*(*(x+r)+j)+t);
                *(*(txtemp)+j) = *(*(*(x+r)+j)+t);
            }

			MatrixMult(xtemp,*K,1,txtemp,*K,xXx);		    /* temp.kk */
			MatrixMult(ytemp,*P,1,txtemp,*K,yXx);			/* temp.pk */
			MatrixMult(ytemp,*P,1,tytemp,*P,yXy);			/* temp.pp */

			for(j=0; j<*K; j++)
			{
				for(jj=0; jj<*K; jj++)
				{
					*(*(M+j)+jj) += *(*(xXx+j)+jj);			/* M */
					*(*(J+j)+jj) += *(*(xXx+j)+jj);			/* J */
				}

				for(i=0; i<*P; i++)
				{
					*(*(G+j)+i) += *(*(xXy+j)+i);			/* G */
				}
			}


			for(i=0; i<*P; i++)
			{
				for(j=0; j<*K; j++)
				{
					*(*(N+i)+j) += *(*(yXx+i)+j);			/* N */
				}
				for(ii=0; ii<*P; ii++)
				{
					*(*(L+i)+ii) += *(*(yXy+i)+ii);			/* L */
					*(*(F+i)+ii) += *(*(yXy+i)+ii);			/* F */
				}
			}

			for(j=0; j<*K; j++)
			{
				*xsingle = *(*(*(x+r)+j)+(t+1));
				for(jj=0; jj<*K; jj++)
				{							                /* H */
					*(*(*(H+j)+jj)) += (*(*(xtemp+jj))) *
						(*xsingle);
				}
				for(i=0; i<*P; i++)
				{							                /* S */
					*(*(*(S+j)+i)) += (*(*(ytemp+i))) *
						(*xsingle);
				}
			}

			for(i=0; i<*P; i++)
			{
				*ysingle = *(*(*(y+r)+i)+t);
				for(j=0; j<*K; j++)					        /* E */
				{
					*(*(*(E+i)+j)) += (*(*(xtemp+j))) *
						(*ysingle);
				}
				*ysingle = *(*(*(y+r)+i)+(t+1));
				for(ii=0; ii<*P; ii++)
				{							                /* Q */
					*(*(*(Q+i)+ii)) += (*(*(ytemp+ii))) *
						(*ysingle);
				}
			}
		}

		/* Add last time T to the J matrix and E list */
  		for(j=0; j<*K; j++)
        {
            *(*(xtemp+j)) = *(*(*(x+r)+j)+(*T-1));
            *(*(txtemp)+j) = *(*(*(x+r)+j)+(*T-1));
        }
        MatrixMult(xtemp,*K,1,txtemp,*K,xXx);
		MatrixSum(J,xXx,J,K,K);
		for(i=0; i<*P; i++)
		{
			for(j=0; j<*K; j++)
			{
				*(*(xtemp+j)) = (*(*(*(x+r)+j)+(*T-1))) * (*(*(*(y+r)+i)+(*T-1)));
			}
			MatrixSum(*(E+i), xtemp, *(E+i), K, one);
		}

        for(j=0; j<*K; j++)
		{
            *(*(M+j)+j) += *(alpha+j);
            *(*(J+j)+j) += *(gamma+j);
        }
        for(i=0; i<*P; i++)
        {
            *(*(L+i)+i) += *(beta+i);
            *(*(F+i)+i) += *(delta+i);
        }
	}

    /* Free Unnecessary Memory */
    /* Keep M, N, L, S, H, J, G, F, E, and Q for later */
	for(i=0; i<*P; i++)
    {
        free(*(ytemp+i));
        free(*(yXx+i));
        free(*(yXy+i));
    }
    for(j=0; j<*K; j++)
    {
        free(*(xtemp+j));
        free(*(xXx+j));
		free(*(xXy+j));
    }
    free(*(txtemp));
    free(*(tytemp));

    free(xsingle);
	free(ysingle);
    free(ytemp);
    free(xtemp);
    free(txtemp);
    free(tytemp);
    free(xXx);
	free(xXy);
    free(yXx);
    free(yXy);

	/* Begin function part 2 (corresponds to simplify in R code) */
	/* Allocate vectors Linv, Minv, Finv, and Jinv */

	Linv = (double**) calloc(*P, sizeof(double*));
	Minv = (double**) calloc(*K, sizeof(double*));
	Finv = (double**) calloc(*P, sizeof(double*));
	Jinv = (double**) calloc(*K, sizeof(double*));
	Nt = (double**) calloc(*K, sizeof(double*));
	Gt = (double**) calloc(*P, sizeof(double*));
	matkk = (double**) calloc(*K, sizeof(double*));
	matpp = (double**) calloc(*P, sizeof(double*));
	matkp = (double**) calloc(*K, sizeof(double*));
	matpk =	(double**) calloc(*P, sizeof(double*));

	for(i=0; i<*P; i++)
	{
		*(Linv+i) = (double*) calloc(*P, sizeof(double));
		*(Finv+i) = (double*) calloc(*P, sizeof(double));
		*(Gt+i) = (double*) calloc(*K, sizeof(double));
		*(matpp+i) = (double*) calloc(*P, sizeof(double));
		*(matpk+i) = (double*) calloc(*K, sizeof(double));
	}
	for(j=0; j<*K; j++)
	{
		*(Minv+j) = (double*) calloc(*K, sizeof(double));
		*(Jinv+j) = (double*) calloc(*K, sizeof(double));
		*(Nt+j) = (double*) calloc(*P, sizeof(double));
		*(matkk+j) = (double*) calloc(*K, sizeof(double));
		*(matkp+j) = (double*) calloc(*P, sizeof(double));
	}

	/* Find matrix inverses */
	MatrixInv(L, *P, Linv, det);							/* L.inv */
	MatrixInv(M, *K, Minv, det);							/* M.inv */
	MatrixInv(F, *P, Finv, det);							/* F.inv */
	MatrixInv(J, *K, Jinv, det);							/* J.inv */

	MatrixTrans(N, Nt, P, K);
	MatrixTrans(G, Gt, K, P);

	MatrixMult(Nt, *K, *P, Linv, *P, matkp);
	MatrixMult(matkp, *K, *P, N, *K, matkk);
	for(j=0; j<*K; j++)
	{
		for(jj=0; jj<*K; jj++)
		{
			*(*(matkk+j)+jj) = *(*(M+j)+jj) - *(*(matkk+j)+jj);
		}
	}
	MatrixInv(matkk, *K, MNLNinv, det);						/* MNLN.inv */

	MatrixMult(N, *P, *K, Minv, *K, matpk);
	MatrixMult(matpk, *P, *K, Nt, *P, matpp);
	for(i=0; i<*P; i++)
	{
		for(ii=0; ii<*P; ii++)
		{
			*(*(matpp+i)+ii) = *(*(L+i)+ii) - *(*(matpp+i)+ii);
		}
	}
	MatrixInv(matpp, *P, LNMNinv, det);						/* LNMN.inv */

	MatrixMult(G, *K, *P, Finv, *P, matkp);
	MatrixMult(matkp, *K, *P, Gt, *K, matkk);
	for(j=0; j<*K; j++)
	{
		for(jj=0; jj<*K; jj++)
		{
			*(*(matkk+j)+jj) = *(*(J+j)+jj) - *(*(matkk+j)+jj);
		}
	}
	MatrixInv(matkk, *K, JGFGinv, det);						/*  JGFG.inv */

	MatrixMult(Gt, *P, *K, Jinv, *K, matpk);
	MatrixMult(matpk, *P, *K, G, *P, matpp);
	for(i=0; i<*P; i++)
	{
		for(ii=0; ii<*P; ii++)
		{
			*(*(matpp+i)+ii) = *(*(F+i)+ii) - *(*(matpp+i)+ii);
		}
	}
	MatrixInv(matpp, *P, FGJGinv, det);						/* FGJG.inv */

	for(j=0; j<*K; j++)
	{
		MatrixMult(Nt, *K, *P, Linv, *P, matkp);
		MatrixMult(matkp, *K, *P, *(S+j), 1, *(HNLS+j));
		MatrixMult(N, *P, *K, Minv, *K, matpk);
		MatrixMult(matpk, *P, *K, *(H+j), 1, *(SNMH+j));
		for(jj=0; jj<*K; jj++)
		{
			*(*(*(HNLS+j)+jj)) = *(*(*(H+j)+jj)) - *(*(*(HNLS+j)+jj));	/* HNLS */
		}
		for(i=0; i<*P; i++)
		{
			*(*(*(SNMH+j)+i)) = *(*(*(S+j)+i)) - *(*(*(SNMH+j)+i));		/* SNMH */
		}
	}

	for(i=0; i<*P; i++)
	{
		MatrixMult(G, *K, *P, Finv, *P, matkp);
		MatrixMult(matkp, *K, *P, *(Q+i), 1, *(EGFQ+i));
		MatrixMult(Gt, *P, *K, Jinv, *K, matpk);
		MatrixMult(matpk, *P, *K, *(E+i), 1, *(QGJE+i));
		for(j=0; j<*K; j++)
		{
			*(*(*(EGFQ+i)+j)) = *(*(*(E+i)+j)) - *(*(*(EGFQ+i)+j));
		}
		for(ii=0; ii<*P; ii++)
		{
			*(*(*(QGJE+i)+ii)) = *(*(*(Q+i)+ii)) - *(*(*(QGJE+i)+ii));
		}
	}

	/* Free Memory */
	for(i=0; i<*P; i++)
	{
		free(*(N+i));
		free(*(L+i));
	    for(j=0; j<*K; j++)
		{
			free(*(*(E+i)+j));
		}
		for(ii=0; ii<*P; ii++)
		{
			free(*(*(Q+i)+ii));
		}
		free(*(E+i));
		free(*(Q+i));
		free(*(F+i));
		free(*(Linv+i));
		free(*(Finv+i));
		free(*(Gt+i));
		free(*(matpp+i));
		free(*(matpk+i));
	}
	free(matpp);
	free(matpk);
	free(Gt);
	free(Linv);
	free(Finv);

	for(j=0; j<*K; j++)
	{
		free(*(M+j));
		for(jj=0; jj<*K; jj++)
		{
			free(*(*(H+j)+jj));
		}
		for(i=0; i<*P; i++)
		{
			free(*(*(S+j)+i));
		}
		free(*(H+j));
		free(*(S+j));
		free(*(J+j));
		free(*(G+j));
		free(*(Minv+j));
		free(*(Jinv+j));
		free(*(Nt+j));
		free(*(matkk+j));
		free(*(matkp+j));
	}
	free(matkk);
	free(matkp);
	free(Nt);
	free(Minv);
	free(Jinv);
	free(E);
    free(Q);
    free(F);
	free(G);
	free(J);
	free(H);
	free(S);
	free(M);
	free(N);
	free(L);
	free(det);
	free(one);
}



/*********************************************************************************/
/*  Simplify code, corresponds to the simplify.2 and simplify functions in R     */
/*  in the case where there are no x's to be estimated.  The inputs are          */
/*  delta = current value of delta, v = current value of d, y, P, T,             */
/*  Rchoice (=a given r if all = FALSE, =R if all = TRUE), all = whether 	     */
/*  estimation is done for a single replicate (FALSE) or all replicates (TRUE),  */
/*  DmeanNox = matrix to hold D mean, DvarNox = matrix to hold Dvar.		     */
/*										                                         */
/*  Note: This is the base variance, to get actual variance for each row, need to*/
/*  multiply by v_i^(-1).							                             */
/*********************************************************************************/

void SimplifyNoX(double *delta, double *v,
	double ***y, int *P, int *T, int *Rchoice, int *all,
	double **DmeanNox, double **DvarNox)
{
	int lower = 0, upper = 0, r, t, i, ii;
	double **p1, **p2, **ytemp, **tytemp, **varp2temp, **yXy, **DmeanNoxt, *det;
    if(*all == 1) {lower = 0; upper = *Rchoice;}
    if(*all == 0) {lower = *Rchoice; upper = (*Rchoice)+1;}

	/* Allocate memory and initialize */
	det = (double*) calloc(1, sizeof(double));
	p1 = (double**) calloc(*P, sizeof(double*));
	p2 = (double**) calloc(*P, sizeof(double*));
	yXy = (double**) calloc(*P, sizeof(double*));
	ytemp = (double**) calloc(*P, sizeof(double*));
	tytemp = (double**) calloc(1, sizeof(double*));
	varp2temp = (double**) calloc(*P, sizeof(double*));
	DmeanNoxt = (double**) calloc(*P, sizeof(double*));
	for(i=0; i<*P; i++)
	{
		*(p1+i) = (double*) calloc(*P, sizeof(double));
		*(p2+i) = (double*) calloc(*P, sizeof(double));
		*(yXy+i) = (double*) calloc(*P, sizeof(double));
		*(ytemp+i) = (double*) calloc(1, sizeof(double));
		*(varp2temp+i) = (double*) calloc(1, sizeof(double));
		*(DmeanNoxt+i) = (double*) calloc(*P, sizeof(double));
		for(ii=0; ii<*P; ii++)
		{
			*(*(p1+i)+ii) = 0;
			*(*(p2+i)+ii) = 0;
		}
	}
	*(tytemp) = (double*) calloc(*P, sizeof(double));

	for(r=lower; r<upper; r++)
	{
		for(t=0; t<((*T)-1); t++)
		{
			for(i=0; i<*P; i++)
			{
				*(*(ytemp+i)) = *(*(*(y+r)+i)+t);
				*(*(tytemp)+i) = *(*(*(y+r)+i)+t);
			}
			MatrixMult(ytemp, *P, 1, tytemp, *P, yXy);

			for(i=0; i<*P; i++)
			{
				for(ii=0; ii<*P; ii++)
				{
					*(*(p1+i)+ii) += *(*(yXy+i)+ii);
				}
			}
			for(i=0; i<*P; i++)
			{
				for(ii=0; ii<*P; ii++)
				{
					*(*(p2+i)+ii) += (*(*(ytemp+i))) *
						(*(*(*(y+r)+ii)+(t+1)));
				}
			}
		}
	}

	/* Begin function, Part 2 */
	/* Corresponds to simplify in R code */

	for(i=0; i<*P; i++)
	{
		*(*(p1+i)+i) += *(delta + i);
	}
	/* This is the base variance, to get actual variance for each row, need to multiply by v_i^(-1) */
	MatrixInv(p1, *P, DvarNox, det);

	/* This is the mean */
	MatrixMult(DvarNox,*P,*P,p2,*P,DmeanNoxt);
	MatrixTrans(DmeanNoxt, DmeanNox, P, P);

	/* Free memory */
	free(*tytemp);
	for(i=0; i<*P; i++)
	{
		free(*(p1+i));
		free(*(p2+i));
		free(*(ytemp+i));
		free(*(yXy+i));
		free(*(varp2temp+i));
		free(*(DmeanNoxt+i));
	}
	free(p1);
	free(p2);
	free(DmeanNoxt);
	free(yXy);
	free(varp2temp);
	free(ytemp);
	free(tytemp);
	free(det);
}





/******************************************************************/
/*  Function to calculate the inverse and determinant of a matrix */
/*  Code thanks to Cherie Ochsenfeld, using LAPACK function       */
/*  "mat" is an n x n matrix, "invmat" is the double pointer to   */
/*  hold the inverse                                              */
/* 								                                  */
/*  Code needs to be compiled with -llapack.                      */
/******************************************************************/

void MatrixInv(double **mat, int n, double **invmat, double *det)
{
	double **u, *w, **v;
	int i, j, k;
	double *fvect, *u2, *vt, *work;
	char jobu = 'A', jobvt = 'A';
	int lwork = -1, info = 0, m = n;

	/* Allocate memory */
	fvect = (double*) calloc(n*n, sizeof(double));
	u2 = (double*) calloc(n*n, sizeof(double));
	w = (double*) calloc(n, sizeof(double));
	vt = (double*) calloc(n*n, sizeof(double));
	work = (double*) calloc(1, sizeof(double));
	u = (double**) calloc(n, sizeof(double*));
	v = (double**) calloc(n, sizeof(double*));

	for(i=0; i<n; i++)
	{
		*(u+i) = (double*) calloc(n, sizeof(double));
		*(v+i) = (double*) calloc(n, sizeof(double));
	}

	/* Convert matrix to vector form for Lapack SVD call */
	for(i=0; i<n; i++) for(j=0; j<n; j++) *(fvect+(j+i*n)) = *(*(mat+j)+i);

	/* Singular value decomposition using LAPACK function */
	F77_CALL(dgesvd)(&jobu,&jobvt,&m,&m,fvect,&m,w,u2,&m,vt,&m,work,&lwork,&info);
	lwork = *work;
	free(work);
	work = (double*) calloc(lwork, sizeof(double));
	F77_CALL(dgesvd)(&jobu,&jobvt,&m,&m,fvect,&m,w,u2,&m,vt,&m,work,&lwork,&info);

	/* Convert U and V to matrix form */
	for(i=0; i<n; i++) for(j=0; j<n; j++) *(*(u+j)+i)=*(u2+(j+i*n));
	for(i=0; i<n; i++) for(j=0; j<n; j++) *(*(v+i)+j)=*(vt+(j+i*n));

	/* Inversion of SVD */
	for(i=0; i<n; i++) for(j=0; j<n; j++) *(*(v+i)+j) = *(*(v+i)+j) *(1/(*(w+j)));

	/* Multiply inverted SVD to get inverted matrix */
	for(i=0; i<n; i++) for(j=0; j<n; j++) for(k=0; k<n; k++) *(*(invmat+i)+j) = *(*(invmat+i)+j) + (*(*(v+i)+k)) * (*(*(u+j)+k));

	/* Calculate the determinate of the matrix */
	*det=0;
	for(i=0; i<n; i++) *det += log(*(w+i));

	/* Free memory */
	for(i=0; i<n; i++)
	{
		free(*(u+i));
		free(*(v+i));
	}
	free(w);
	free(u);
	free(v);
	free(fvect);
	free(u2);
	free(vt);
	free(work);
}

/****************************************************************************/
/*  This is a function to perform matrix multiplication, m1 %*% m2 = sol.   */
/*  This code is thanks to Cherie Ochsenfeld, and uses the BLAS function.   */
/*  m1r and m1c are the number of rows and cols of m1, and m2c is the       */
/*  number of columns of m2. 					`	    */
/****************************************************************************/

void MatrixMult(double **m1, int m1r, int m1c, double **m2,
	int m2c, double **sol)
{
	int i,j;
	double *A, *B, *C, alph=1.0, bta=0.0;
	char transa='N', transb='N';

	/* Allocate memory */
	A = (double*) calloc((m1r*m1c), sizeof(double));
	B = (double*) calloc((m1c*m2c), sizeof(double));
	C = (double*) calloc((m1r*m2c), sizeof(double));

	/* Turn matrix into a vector */
	for(i=0; i<m1c; i++) for(j=0; j<m1r; j++) *(A+(i*m1r+j))=*(*(m1+j)+i);
	for(i=0; i<m2c; i++) for(j=0; j<m1c; j++) *(B+(i*m1c+j))=*(*(m2+j)+i);

	/* Call BLAS function */
	F77_CALL(dgemm)(&transa, &transb, &m1r, &m2c, &m1c, &alph, A, &m1r, B, &m1c, &bta, C, &m1r);

	/* set the solution as a matrix */
	for(i=0; i<m2c; i++) for(j=0; j<m1r; j++) *(*(sol+j)+i) = *(C+(i*m1r+j));

	/* Free memory */
	free(A);
	free(B);
	free(C);
}


/************************************************************************/
/*  Function to perform element-wise addition of two matrices,          */
/*  m1 + m2 = sum.                                                      */
/************************************************************************/

void MatrixSum(double **m1, double **m2, double **sum, int *row, int *col)
{
	int i, j;
	for(i=0; i<*row; i++)
	{
		for(j=0; j<*col; j++)
		{
			*(*(sum+i)+j) = (*(*(m1+i)+j)) + (*(*(m2+i)+j));
		}
	}
}

/*************************************************************/
/*  Function to transpose a matrix contained in a double     */
/*  pointer, mt = t(m).  "row" and "col" are the number of   */
/*  rows and columns of the original matrix, m.              */
/*************************************************************/

void MatrixTrans(double **m, double **mt, int *row, int *col)
{
	int i,j;
	for(i=0; i<*row; i++)
	{
		for(j=0; j<*col; j++)
		{
			*(*(mt+j)+i) = *(*(m+i)+j);
		}
	}
}
