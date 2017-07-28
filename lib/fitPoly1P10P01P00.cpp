//
//  fitPoly.cpp
//  fits polynomial of degree 2 to data
//
//  Created by Friedhelm Serwane on 5/21/13.
//  Copyright (c) 2013 Friedhelm Serwane. All rights reserved.
//

//

#include <iostream>
#include <sstream>
#include<fstream>
#include<cmath>
#include<vector>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <matrix.h>
#include <mex.h>




using namespace std;

struct data {
    size_t n;
    double * z;
    double * SX;
    double * SY;
    double * weights;
};

int modelFunction_f (const gsl_vector * x, void *data, gsl_vector * f)
{ size_t n = ((struct data *)data)-> n;
    double *SX = ((struct data *)data)->SX;
    double *SY = ((struct data *)data)->SY;
    double *z = ((struct data *)data)->z;
    double *weights = ((struct data *)data)->weights;
    
  
    double p10=gsl_vector_get(x,0);
    double p01=gsl_vector_get(x, 1);
    double p00=gsl_vector_get(x, 2);
    
    size_t i;
    for (i = 0; i < n; i++)
    {
        /* Model 2 Yi = p20 x^2 +p02 y^2 +p11 x y + p10 x +p01 y + p00 */
        double tx = SX[i];
        double ty = SY[i];
        double Zi = p10*tx+p01*ty+p00;
        gsl_vector_set (f, i, (Zi-z[i])/weights[i]);
    }
    return GSL_SUCCESS;
}

int
modelFunction_df (const gsl_vector * x, void *data, gsl_matrix * J)
{
    size_t n = ((struct data *)data)->n;
    double *SX = ((struct data *)data)->SX;
    double *SY = ((struct data *)data)->SY;
    double *weights = ((struct data *)data)->weights;

    
    double p10=gsl_vector_get(x, 0);
    double p01=gsl_vector_get(x, 1);
    double p00=gsl_vector_get(x, 2);
   
    size_t i;
    for (i = 0; i < n; i++)
    {
        /* Jacobian matrix J(i,j) = dfi / dxj, */
        /* Zi = p11 x y + p10 x +p01 y + p00*/
        /* and the xj are the parameters (p20, p02,p11,p10,p01,p00) */
        double tx = SX[i];
        double ty = SY[i];
        double tw = weights[i];
        
       
        
        gsl_matrix_set (J, i, 0, 1/tw*(tx));
        gsl_matrix_set (J, i, 1, 1/tw*(ty));
        gsl_matrix_set (J, i, 2, 1/tw*(1));
    
    }
    return GSL_SUCCESS;
}

int
modelFunction_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
    modelFunction_f (x, data, f);
    modelFunction_df (x, data, J);
    return GSL_SUCCESS;
}

//The main part of the program sets up a Levenberg-Marquardt solver. 
//void print_state (size_t iter, gsl_multifit_fdfsolver * s);

double *fPoly (double fitParInitial[],double dataX[],double dataY[], double dataZ[],double weights[],int dimData,int dimfpar)
{
   
    unsigned int i=0, iter = 0;
    string line, temp;
  
    const size_t p = 3;
    //p= dimfpar;
    int lout=2*p;
    #define N dimData
    //#define ARRAY_SIZE(array) (sizeof((array))/sizeof((array[0])))
   
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    int status;
    const size_t n = N;
    
    int j, k;
    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    
    double* z = new double[N];
    double* SX = new double[N];
    double* SY = new double[N];
    double* SW = new double[N];
    
    
   
    struct data d = { n, z, SX,SY,SW};
    gsl_multifit_function_fdf f;
    double x_init[p] = {0.0};
    
    for (i=0; i<p; i++) {
        x_init[i]=fitParInitial[i];
    }
    gsl_vector_view x = gsl_vector_view_array (x_init, p);
    f.f = &modelFunction_f;
    f.df = &modelFunction_df;
    f.fdf = &modelFunction_fdf;
    f.n = n;
    f.p = p;
    f.params = &d;
    
    for (i=0; i<N; i++) {
        z[i]=dataZ[i];
        SX[i]=dataX[i];
        SY[i]=dataY[i];
        SW[i]=weights[i];
    }
    
 
    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc (T, n, p);
    gsl_multifit_fdfsolver_set (s, &f, &x.vector);
    //print_state (iter, s);
    do
    {
        iter++;
        status = gsl_multifit_fdfsolver_iterate (s);
        //printf ("status = %s\n", gsl_strerror (status));
        //print_state (iter, s);
        if (status)
            break;
        status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
        //cout<<iter;
    }
    while (status == GSL_CONTINUE && iter < 10000);
    gsl_multifit_covar (s->J, 0.0, covar);
    #define FIT(i) gsl_vector_get(s->x, i)
    #define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof));
    double rsquare=pow(chi, 2.0) / dof;
        //printf("chisq/dof = %g\n", pow(chi, 2.0) / dof);
        //printf ("sigma = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
        //printf ("mu = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
    // printf ("status = %s\n", gsl_strerror (status));
    
    int lfitPar=p;
    //generate results: 0-5 fitparameters 5
    
    static double fitResults[7]={0.0};
    iter=0;
    do {
        fitResults[iter]=gsl_vector_get(s->x, iter);
        fitResults[iter+lfitPar]=c*ERR(iter);
        iter++;
    } while (iter<lfitPar);
    
   
    fitResults[2*lfitPar]=rsquare;

    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);
    /*cout << rsquare;*/
    return fitResults;
}

void
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
    printf ("iter: %3u x = % 15.8f % 15.8f "
            "|f(x)| = %g\n",
            iter,
            gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1),
            gsl_blas_dnrm2 (s->f));
    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //declare variables
    mxArray *fparIn_in_m, *dataX_in_m, *dataY_in_m,*dataZ_in_m, *fpar_out_m, *errfpar_out_m,*chisq_out_m, *weights_in_m;
    const mwSize *dimsfparIn;
    const mwSize *dimsOut;
    const mwSize *dimsDataX;
    const mwSize *dimsDataY;
    const mwSize *dimsDataZ;
    const mwSize *dimsWeights;
    
    double *fparIn, *fpar, *errfpar, *chisq, *dataX,*dataY,*dataZ,*weights;
    
    int dimfparIn, dimDataX,dimDataY,dimDataZ,dimWeights;
    int numdimsfparIn,numdimsDataX,numdimsDataY,numdimsDataZ,numdimsWeights;
    int numdimsOut, dimxOut,dimyOut;
    int i,j;
    
    
    //associate inputs
    fparIn_in_m = mxDuplicateArray(prhs[0]);
    dataX_in_m = mxDuplicateArray(prhs[1]);
    dataY_in_m = mxDuplicateArray(prhs[2]);
    dataZ_in_m = mxDuplicateArray(prhs[3]);
    weights_in_m = mxDuplicateArray(prhs[4]);
    
    //figure out dimensions
    //inputs
    dimsfparIn = mxGetDimensions(prhs[0]);
    numdimsfparIn = mxGetNumberOfDimensions(prhs[0]);
    dimfparIn = (int)dimsfparIn[0]; 
    int lfpar=3;
        
    dimsDataX=mxGetDimensions(prhs[1]);
    numdimsDataX = mxGetNumberOfDimensions(prhs[1]);
    dimDataX = (int)dimsDataX[0]; 
    
    
    dimsDataY=mxGetDimensions(prhs[2]);
    numdimsDataY = mxGetNumberOfDimensions(prhs[2]);
    dimDataY = (int)dimsDataY[0];
    
    dimsDataZ=mxGetDimensions(prhs[3]);
    numdimsDataZ = mxGetNumberOfDimensions(prhs[3]);
    dimDataZ = (int)dimsDataZ[0];
    
    dimsWeights=mxGetDimensions(prhs[4]);
    numdimsWeights = mxGetNumberOfDimensions(prhs[4]);
    dimWeights = (int)dimsWeights[0];

    

    //associate outputs
    fpar_out_m = plhs[0] = mxCreateDoubleMatrix(1,lfpar,mxREAL);
    errfpar_out_m = plhs[1] = mxCreateDoubleMatrix(1,lfpar,mxREAL);
    chisq_out_m = plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    //associate pointers
    fparIn = mxGetPr(fparIn_in_m);
    dataX= mxGetPr(dataX_in_m);
    dataY= mxGetPr(dataY_in_m);
    dataZ= mxGetPr(dataZ_in_m);
    weights= mxGetPr(weights_in_m);
    
    fpar= mxGetPr(fpar_out_m);
    errfpar= mxGetPr(errfpar_out_m);
    chisq=mxGetPr(chisq_out_m);
    //do something
    //for(i=0;i<dimx;i++)
    //{
    //    for(j=0;j<dimy;j++)
    //    {
    //        //mexPrintf("element[%d][%d] = %f\n",j,i,a[i*dimy+j]);
            
    //        d[i*dimy+j] = a[i*dimy+j]*a[i*dimy+j]; //squares a
            
    //    }
    //}
    //int pfval=fGauss();
    //int fval=*pfval;
    
    
    double *p1;
 
    p1=fPoly(fparIn,dataX,dataY,dataZ,weights,dimDataX,dimfparIn);
        for (i=0;i < lfpar;i++)
    {
        fpar[i]=*(p1+i);
        errfpar[i]=*(p1+i+lfpar);
    }
    chisq[0]=*(p1+2*lfpar);
    //r[0]=4;
    //mexPrintf("fit results");
    return;
}

