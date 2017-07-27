//  fitGauss.cpp
//  fitGauss
//
//  Created by Friedhelm Serwane on 5/21/13.
//  Copyright (c) 2013 Friedhelm Serwane. All rights reserved.
//  based on this example: http://www.gnu.org/software/gsl/manual/gsl-ref_38.html#SEC518



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

//order of fit paramters: sigma mu amp offset


using namespace std;


struct data {
    size_t n;
    double * y;
    double * SB;
};

int expb_f (const gsl_vector * x, void *data, gsl_vector * f)
{ size_t n = ((struct data *)data)-> n;
  double *SB = ((struct data *)data)->SB;
  double *y = ((struct data *)data)->y;
  
  double sigma = gsl_vector_get(x,0);
  double mu = gsl_vector_get(x, 1);
  double amplitude=gsl_vector_get(x, 2);
  double offset=gsl_vector_get(x, 3);
  
  size_t i;
  for (i = 0; i < n; i++)
  {
      //Here comes the fit model
      /* Model Yi = (1/sqrt(2*M_PI*sigma*sigma))*exp(-(t-mu)*(t-mu)/(2*sigma*sigma)) */
      /* Model 2 Yi = amplitude*exp(-(t-mu)*(t-mu)/(2*sigma*sigma))+offset */
      double t = SB[i];
      double Yi = amplitude*exp(-(t-mu)*(t-mu)/(2*sigma*sigma))+offset ;
      gsl_vector_set (f, i, Yi-y[i]);
  }
  return GSL_SUCCESS;
}

int
        expb_df (const gsl_vector * x, void *data, gsl_matrix * J)
{
    size_t n = ((struct data *)data)->n;
    double *SB = ((struct data *)data)->SB;
    
    double sigma = gsl_vector_get (x, 0);
    double mu = gsl_vector_get (x, 1);
    double amplitude=gsl_vector_get(x, 2);
    double offset=gsl_vector_get(x, 3);
    
    size_t i;
    for (i = 0; i < n; i++)
    {
        /* Jacobian matrix J(i,j) = dfi / dxj, */
        /* Yi = amplitude*exp(-(t-mu)*(t-mu)/(2*sigma*sigma))+offset*/
        /* and the xj are the parameters (sigma, mu) */
        double t = SB[i];
        double e = exp(-(t-mu)*(t-mu)/(2*sigma*sigma));
        
        gsl_matrix_set (J, i, 0, (amplitude*e*((t-mu)*(t-mu)/(sigma*sigma*sigma))));
        
        gsl_matrix_set (J, i, 1, ((amplitude * e*((t-mu)/(sigma*sigma)) )));
        gsl_matrix_set (J, i, 2, (e ));
        gsl_matrix_set (J, i, 3, (1));
        
        /*gsl_matrix_set (J, i, 1, (-(A/(sigma*sigma) * e * (t-mu))));*/
    }
    return GSL_SUCCESS;
}

int
        expb_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)
{
    expb_f (x, data, f);
    expb_df (x, data, J);
    return GSL_SUCCESS;
}

//The main part of the program sets up a Levenberg-Marquardt solver.
void print_state (size_t iter, gsl_multifit_fdfsolver * s);


double *fGauss (double fitParInitial[],double dataX[],double dataY[],int dimData)
{
    unsigned int i=0, iter = 0;
    string line, temp;
    
    //ifstream read;
    //read.open ("fake_Gaus_data.txt");
    // if(read.is_open())
    //     while(!read.eof()){
    //         getline(read,line);
    //         i++;}
    // read.close();
    //   cout << sizeof(dataX);
    //   cout << sizeof(dataY);
    #define N dimData
    //#define ARRAY_SIZE(array) (sizeof((array))/sizeof((array[0])))
    // int N=ARRAY_SIZE(dataX);
    
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    int status;
    const size_t n = N;
    const size_t p = 4;
    int j, k;
    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    double y[N], SB[N];
    struct data d = { n, y, SB};
    gsl_multifit_function_fdf f;
    double x_init[4] = {0.0};
    
    for (i=0; i<p; i++) {
        x_init[i]=fitParInitial[i];
    }
    
    
    gsl_vector_view x = gsl_vector_view_array (x_init, p);
    
    
    
    f.f = &expb_f;
    f.df = &expb_df;
    f.fdf = &expb_fdf;
    f.n = n;
    f.p = p;
    f.params = &d;
    
    
    for (i=0; i<N; i++) {
        y[i]=dataY[i];
        SB[i]=dataX[i];
    }
    
    
    //for (i = 0; i < n; i++)
    //{
    //    printf ("data: %g %g\n", SB[i], y[i]);
    //}
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
    }
    while (status == GSL_CONTINUE && iter < 5000);
    gsl_multifit_covar (s->J, 0.0, covar);
    #define FIT(i) gsl_vector_get(s->x, i)
    #define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c =chi / sqrt(dof);
    double rsquare=pow(chi, 2.0) / dof;
    //printf("chisq/dof = %g\n", pow(chi, 2.0) / dof);
    //printf ("sigma = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
    //printf ("mu = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
    
    // printf ("status = %s\n", gsl_strerror (status));
    
    
    
    int lfitPar=p;
    //generate results: 0-4 fitparameters 5
    static double fitResults[9]={0.0};
    iter=0;
    do {
        fitResults[iter]=gsl_vector_get(s->x, iter);
        fitResults[iter+lfitPar]=c*ERR(iter);
        iter++;
    } while (iter<lfitPar);
    
    
    
    fitResults[8]=rsquare;
    
    
    
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
    mxArray *a_in_m, *dataX_in_m, *dataY_in_m, *fpar_out_m, *errfpar_out_m,*chisq_out_m;
    const mwSize *dims;
    const mwSize *dimsOut;
    const mwSize *dimsDataX;
    const mwSize *dimsDataY;
    double *a, *fpar, *errfpar, *chisq, *dataX,*dataY;
    int dimx, dimy,dimxDataX,dimyDataX,dimxDataY,dimyDataY, numdims,numdimsDataX,numdimsDataY;
    int numdimsOut, dimxOut,dimyOut;
    int i,j;
    
    //associate inputs
    a_in_m = mxDuplicateArray(prhs[0]);
    dataX_in_m = mxDuplicateArray(prhs[1]);
    dataY_in_m = mxDuplicateArray(prhs[2]);
    
    //figure out dimensions
    //inputs
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dimy = (int)dims[0]; dimx = (int)dims[1];
    
    
    dimsDataX=mxGetDimensions(prhs[1]);
    numdimsDataX = mxGetNumberOfDimensions(prhs[1]);
    dimyDataX = (int)dimsDataX[0]; dimxDataX = (int)dimsDataX[1];
    
    
    dimsDataY=mxGetDimensions(prhs[2]);
    numdimsDataY = mxGetNumberOfDimensions(prhs[2]);
    dimyDataY = (int)dimsDataY[0]; dimxDataY = (int)dimsDataY[1];
    
    //Outputs
    //dimsOut = mxGetDimensions(plhs[0]);
    //numdimsOut = mxGetNumberOfDimensions(plhs[0]);
    //dimyOut = (int)dimsOut[0]; dimxOut = (int)dimsOut[1];
    
    
    //associate outputs
    fpar_out_m = plhs[0] = mxCreateDoubleMatrix(1,4,mxREAL);
    errfpar_out_m = plhs[1] = mxCreateDoubleMatrix(1,4,mxREAL);
    chisq_out_m = plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    //associate pointers
    a = mxGetPr(a_in_m);
    dataX= mxGetPr(dataX_in_m);
    dataY= mxGetPr(dataY_in_m);
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
    p1=fGauss(a,dataX,dataY,dimxDataX);
    int nofFpar=4;
    for (i=0;i < nofFpar;i++)
    {
        fpar[i]=*(p1+i);
        errfpar[i]=*(p1+i+nofFpar);
    }
    chisq[0]=*(p1+8);
    //r[0]=4;
    //mexPrintf("fit results");
    return;
}

