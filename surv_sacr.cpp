//
//  The primal-dual interior point algorithm for the survival-sacrifice model
//
//  Created by Piet Groeneboom on 10/03/18.
//  Copyright (c) 2018 Piet Groeneboom. All rights reserved.
//

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string.h>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

#define SQR(x) ((x)*(x))

double  *d,*DD,**L,**b;
double	*grad,*g,*labda,*w,*D,*delta_z,*delta_labda,*delta_w,*GtLabda;
double	*z1,*labda1,*w1,*g1,*Gtranspose_D_gz,*Gtranspose_Winverse_E,*G_delta_z;
double	*A,**Hessian,**Gtranspose_D_G;

typedef struct
{
    double x;
    int delta1;
    int delta2;
}
data_object;

void    Cholesky_dec(int n, double **b, double D[], double **L);
void    Cholesky_sol(int n, double **b, double z[], double x[], double *D, double **L);
void    sort_ties(int n, int delta1[], int delta2[], double xx[],
               double xx_new[], int **freq, int *n_new);
void    sort_data(int n, double xx[], int delta1[], int delta2[]);
int     CompareTime(const void *a, const void *b);

double  LoopInterior(int m, int n, int n2, double mu, double sigma, double beta, int freq[],double z[]);
void	Compute_z(int n1, double p[], double z[]);
void 	InteriorPoint(int n, double xcoor[], double z[], int **freq);
void	Compute_delta_z(int n, int n2, double **Hessian, double A[], double delta_z[]);
void 	Compute_gradient(int n, int **freq, double x[], double grad[]);
void 	initializeInterior(int n, double z[]);
void 	Compute_Gtranspose_Labda(int n, int n2, double labda[], double GtLabda[]);
void 	Compute_Gtranspose_D_gz(int n, int n2, double d[], double g[], double Gtranspose_D_gz[]);
void 	Compute_Gtranspose_Winverse_E(int n, int n2, double w[], double Gtranspose_Winverse_E[]);
void 	Compute_g(int n, int n2, double z[], double g[]);
void	Compute_G_delta_z(int n, int n2, double delta_z[], double G_delta_z[]);
void	Compute_A(int m, double mu, double sigma, double z[], double w[], double Gtranspose_D_gz[],
                  double grad[], double Gtranspose_Winverse_E[], double A[]);
void	Compute_delta_labda(int m, double mu, double sigma, double G_delta_z[], double g[],
                            double d[], double w[], double delta_labda[]);
void	Compute_delta_w(int m, double mu, double sigma, double labda[], double w[], double delta_labda[],
                        double D[], double delta_w[]);
double 	phiInterior(int n, int n2, int **freq, double z[]);
void	gradientInterior(int n, int **freq, double z[], double GtLabda[], double grad[]);
void	Compute_Hessian(int n, int **freq, double z[], double **Gtranspose_D_G, double **Hessian);
void	Compute_D(int m, double labda[], double w[], double d[]);
void	Compute_Gtranspose_D_G(int n, int n2, double d[], double **Gtranspose_D_G);
double  dotest(int n, int m, double z[], int **freq, double labda[], double w[]);
double 	ComputeFeasibleStep(int n, int n2, int m, double z1[], int **freq,
                            double delta_z[], double labda1[], double delta_labda[], double w1[], double delta_w[]);
void	transfer(int n, double a[], double b[]);




// [[Rcpp::export]]

List Compute_estimates(DataFrame input)
{
    int			ndata,m1,m2,**freq;
	int 		  i,j,n,n2,*delta1,*delta2,*index1,*index2;
	double		*xcoor,*mle,*mle1,*mle2,*xx;
    
    DataFrame DF = Rcpp::DataFrame(input);
    NumericVector xcoor0 = DF["V1"];
    IntegerVector delta01 = DF["V2"];
    IntegerVector delta02 = DF["V3"];
    
    // determine the sample size
    
    ndata = (int)xcoor0.size();
    
    delta1 = new int[ndata+1];
    delta2 = new int[ndata+1];
    
    xcoor = new double[ndata+2];
    
    for (i=1;i<=ndata;i++)
    {
        xcoor[i]=(double)xcoor0[i-1];
        delta1[i]=(int)delta01[i-1];
        delta2[i]=(int)delta02[i-1];
    }
 
    n2=2*ndata;
    
    freq = new int *[ndata+1];
    for (i=0;i<ndata+1;i++)
        freq[i] = new int[4];
    
    xx = new double[ndata+2];
    
    //the following matrices are used in theh Cholesky decomposition
    
    L = new double *[n2+1];
    for (i=0;i<n2+1;i++)
        L[i] = new double [3];
    
    for (i=0;i<=n2;i++)
    {
        L[i][0]=1;
        L[i][1]=L[i][2]=0;
    }
    
    d= new double[n2+1];
    DD= new double[n2+1];
    
    
    b = new double *[n2+1];
    for (i=0;i<n2+1;i++)
        b[i] = new double [3];
    
    for (i=0;i<n2+1;i++)
        d[i]=0;
    
    
    mle = new double[n2+1];
    
    xcoor[0]=0;
    
    // the "point at infinity":
    xcoor[ndata+1]=10000;
    
    // put the data in order w.r.t. the ordering of xcoor
    
    sort_data(ndata,xcoor,delta1,delta2);
    
    //  sort data with ties; the observation points xcoor[i] are replaced
    //  by the strictly different points xx[i]
    
    sort_ties(ndata,delta1,delta2,xcoor,xx,freq,&n);
    
    // the number of strictly different observations, seen if following statement is not commented out
    
    // printf("\nn = %5d\n\n",n);
    
    // Do the interior point method
    
    InteriorPoint(n,xx,mle,freq);
    
    // create right-continuous functions on all points from the solution
    
    index1 = new int[n+1];
    index2 = new int[n+1];
    mle1 = new double[n+1];
    mle2 = new double[n+1];
    
    
    j=0;
    
    for (i=1;i<=n;i++)
    {
        if (freq[i][3]>0 || freq[i][2]>0)
        {
            j++;
            index1[j]=i;
        }
    }
    
    m1=j;
    
    for (j=0;j<index1[1];j++)
        mle1[j]=0;
    
    for (i=1;i<m1;i++)
    {
        for (j=index1[i];j<index1[i+1];j++)
            mle1[j]=mle[2*index1[i]-1];
    }
    
    for (j=index1[m1];j<=n;j++)
        mle1[j]=mle[2*index1[m1]-1];
    
    j=0;
    
    for (i=1;i<=n;i++)
    {
        if (freq[i][1]>0 || freq[i][2]>0)
        {
            j++;
            index2[j]=i;
        }
    }
    
    m2=j;
    
    for (j=0;j<index2[1];j++)
        mle2[j]=0;
    
    for (i=1;i<m2;i++)
    {
        for (j=index2[i];j<index2[i+1];j++)
        {
            mle2[j]=mle[2*index2[i]];
            if (mle2[j]<mle1[j])
                mle2[j]=mle1[j];
        }
    }
    
    for (j=index2[m2];j<=n;j++)
    {
        mle2[j]=mle[2*index2[m2]];
        if (mle2[j]<mle1[j])
            mle2[j]=mle1[j];
    }
    
    NumericMatrix out0 = NumericMatrix(n,2);
    
    for (i=0;i<n;i++)
    {
        out0(i,0)=xx[i+1];
        out0(i,1)=mle1[i+1];
    }
    
    NumericMatrix out1 = NumericMatrix(n,2);
    
    for (i=0;i<n;i++)
    {
        out1(i,0)=xx[i+1];
        out1(i,1)=mle2[i+1];
    }
    
    // make the list for the output, containing the two estimates
    
    List out = List::create(Rcpp::Named("MLE1")=out0,Rcpp::Named("MLE2")=out1);

    // free memory

    delete[] delta1, delete[] delta2, delete[] xcoor, delete[] xx, delete[] mle;
    
    delete[] mle1, delete[] mle2, delete[] index1, delete[] index2;
    
    for (i=0;i<n2+1;i++)
        delete[] L[i];
    delete[] L;
    
    for (i=0;i<n2+1;i++)
        delete[] b[i];
    delete[] b;
    
    for (i=0;i<ndata+1;i++)
        delete[] freq[i];
    delete[] freq;
    
    delete[] d; delete[] DD;
    
    return out;
    
}

void Cholesky_dec(int n, double **b, double DD[], double **L)
{
    int i;
    
    // b is an a version of the matrix B in (2.29) on p. 25 of Green and Silverman (1994)
    
    DD[1]= b[1][0];
    L[2][1] = b[2][1]/DD[1];
    DD[2]= b[2][0]-SQR(L[2][1])*DD[1];
    
    for (i=3;i<=n;i++)
    {
        L[i][2] = b[i][2]/DD[i-2];
        L[i][1] = (b[i][1]-L[i-1][1]*L[i][2]*DD[i-2])/DD[i-1];
        DD[i]= b[i][0]-SQR(L[i][1])*DD[i-1]-SQR(L[i][2])*DD[i-2];
    }
    
}

void Cholesky_sol(int n, double **b, double z[], double x[], double *DD, double **L)
{
    int i;
    double *u,*v;
    
    
    u= new double[n+1];
    v= new double[n+1];
    
    Cholesky_dec(n,b,DD,L);
    
    u[1]=z[1];
    u[2]=z[2]-L[2][1]*u[1];
    
    for (i=3;i<=n;i++)
        u[i]=z[i]-L[i][1]*u[i-1]-L[i][2]*u[i-2];
    
    for (i=1;i<=n;i++)
        v[i]=u[i]/DD[i];
    
    x[n]=v[n];
    x[n-1]=v[n-1]-L[n][1]*x[n];
    
    for (i=n-2;i>=1;i--)
        x[i]=v[i]-L[i+1][1]*x[i+1]-L[i+2][2]*x[i+2];
    
    delete[] u;
    delete[] v;
    
}

void sort_ties(int n, int delta1[], int delta2[], double xx[],
               double xx_new[], int **freq, int *n_new)
{
    int i,j;
    
    for (i=1;i<=n;i++)
        for (j=1;j<=3;j++)
            freq[i][j]=0;
    
    j=1;
    
    if (delta1[1]==0 && delta2[1]==0)
        freq[1][1]++;
    else
    {
        if (delta1[1]==1 && delta2[1]==0)
            freq[1][2]++;
        else
            freq[1][3]++;
    }
    
    xx_new[1]=xx[1];
    
    for (i=2;i<=n;i++)
    {
        if (xx[i]>xx[i-1])
        {
            j++;
            xx_new[j]=xx[i];
        }
        
        if ((delta1[i]==0) && (delta2[i]==0))
            freq[j][1]++;
        else
        {
            if (delta1[i]==1 && delta2[i]==0)
                freq[j][2]++;
            else
                freq[j][3]++;
        }
    }
    *n_new=j;
}

int CompareTime(const void *a, const void *b)
{
    if ((*(data_object *) a).x < (*(data_object *) b).x)
        return -1;
    if ((*(data_object *) a).x > (*(data_object *) b).x)
        return 1;
    return 0;
}

void sort_data(int n, double xx[], int delta1[], int delta2[])
{
    int i;
    data_object *obs;
    
    obs= new data_object[n];
    
    for (i=0;i<n;i++)
    {
        obs[i].x=xx[i+1];
        obs[i].delta1=delta1[i+1];
        obs[i].delta2=delta2[i+1];
    }
    
    qsort(obs,n,sizeof(data_object),CompareTime);
    
    for (i=0;i<n;i++)
    {
        xx[i+1]=obs[i].x;
        delta1[i+1]=obs[i].delta1;
        delta2[i+1]=obs[i].delta2;
    }
    
    delete[] obs;
}


double	DoNeighborhoodTest(int m, int n, int n2, int **freq, double mu, double beta,double z[])
{
    int i;
    double ladba_times_w,sum=0.0,sum1=0.0,sum2=0.0,mu1,norm,gamma=0.01,chi=0.7;
    
    GtLabda= new double[n2+1];
    
    Compute_Gtranspose_Labda(n,n2,labda1,GtLabda);
    gradientInterior(n,freq,z1,GtLabda,grad);
    
    Compute_g(n,n2,z1,g1);
    
    sum1=0;
    for (i=1;i<=n2;i++)
        sum1=sum1+SQR(grad[i]);
    sum1=sqrt(sum1);
    
    sum2=0;
    for (i=1;i<=m;i++)
        sum2=sum2+SQR(g1[i]+w1[i]);
    sum2=sqrt(sum2);
    
    norm=fmax(sum1,sum2);
    
    ladba_times_w=labda1[1]*w1[1];
    sum=ladba_times_w;
    for (i=2;i<=m;i++)
    {
        sum=sum+labda1[i]*w1[i];
        if (ladba_times_w>labda1[i]*w1[i])
            ladba_times_w=labda1[i]*w1[i];
    }
    
    mu1=sum/m;
    
    while (((norm>beta*mu) || (ladba_times_w<gamma*mu) || (mu1>(1-0.01*chi)*mu)) && (chi>1.0e-1))
    {
        chi = 0.8*chi;
        
        for (i=1;i<=n2;i++)
            z1[i]=z[i]+ chi*delta_z[i];
        
        for (i=1;i<=m;i++)
        {
            labda1[i]=labda[i]+ chi*delta_labda[i];
            w1[i]=w[i]+ chi*delta_w[i];
        }
        
        
        Compute_Gtranspose_Labda(n,n2,labda1,GtLabda);
        gradientInterior(n,freq,z1,GtLabda,grad);
        Compute_g(n,n2,z1,g1);
        
        sum1=0;
        for (i=1;i<=n2;i++)
            sum1=sum1+SQR(grad[i]);
        sum1=sqrt(sum1);
        
        sum2=0;
        for (i=1;i<=m;i++)
            sum2=sum2+SQR(g1[i]+w1[i]);
        sum2=sqrt(sum2);
        
        norm=fmax(sum1,sum2);
        
        sum=labda1[1]*w1[1];
        for (i=2;i<=m;i++)
        {
            sum=sum+labda1[i]*w[i];
            if (ladba_times_w>labda1[i]*w1[i])
                ladba_times_w=labda1[i]*w1[i];
        }
        
        mu1=sum/m;
    }
    
    transfer(n2,z1,z);
    transfer(m,labda1,labda);
    transfer(m,w1,w);
    
    return mu1;
}

double LoopInterior(int m, int n, int n2, double mu, double sigma, double beta, int **freq,double z[])
{
    double a,mu1;
    
    Compute_g(n,n2,z,g);
    
    Compute_D(m,labda,w,D);
    
    Compute_Gtranspose_Labda(n,n2,labda,GtLabda);
    gradientInterior(n,freq,z,GtLabda,grad);
    
    Compute_Gtranspose_D_G(n,n2,D,Gtranspose_D_G);
    Compute_Hessian(n,freq,z,Gtranspose_D_G,Hessian);
    
    Compute_Gtranspose_Winverse_E(n,n2,w,Gtranspose_Winverse_E);
    Compute_Gtranspose_D_gz(n,n2,D,g,Gtranspose_D_gz);
    Compute_A(n2,mu,sigma,z,w,Gtranspose_D_gz,grad,Gtranspose_Winverse_E,A);
    
    // Computes expression on the right of first equation in (7.21), called A
    
    Compute_delta_z(n,n2,Hessian,A,delta_z);
    
    //Cholesky_sol(n,Hessian,A,delta_z,D,L);
    
    Compute_G_delta_z(n,n2,delta_z,G_delta_z);
    Compute_delta_labda(m,mu,sigma,G_delta_z,g,D,w,delta_labda);
    Compute_delta_w(m,mu,sigma,labda,w,delta_labda,D,delta_w);
    
    transfer(n2,z,z1);
    transfer(m,labda,labda1);
    transfer(m,w,w1);
    
    a=ComputeFeasibleStep(n,n2,m,z1,freq,delta_z,labda1,delta_labda,w1,delta_w);
    
    mu1 = DoNeighborhoodTest(m,n,n2,freq,mu,beta,z);
    return mu1;
}


void 	InteriorPoint(int n, double xx[], double z[], int **freq)
{
    double			a,mu=0.5;
    double			phinew,sum1,sum2,norm;
    int 			m,n2,i,j;
    double			tol=1.0e-10;
    double			beta,sigma=0.5,sum;
    
    
    
    beta=10e7;
    norm=10;
    sum=10;
    
    m=3*n;
    n2=2*n;
    
    z1= new double[2*n+2];
    g = new double[m+1];
    g1 = new double[m+1];
    Gtranspose_D_gz = new double[n2+1];
    Gtranspose_Winverse_E = new double[n2+1];
    G_delta_z = new double[m+1];
    
    grad = new double[n2+1];
    A = new double[n2+1];
    D = new double[m+1];
    labda = new double[m+1];
    w = new double[m+1];
    labda1 = new double[m+1];
    w1 = new double[m+1];
    delta_z = new double[m+1];
    delta_labda = new double[m+1];
    delta_w = new double[m+1];
    GtLabda = new double[n2+1];
    
    Hessian = new double *[n2+1];
    for (i=0;i<n2+1;i++)
        Hessian[i] = new double [3];
    
    Gtranspose_D_G = new double *[n2+1];
    for (i=0;i<n2+1;i++)
        Gtranspose_D_G[i] = new double [3];
    
    initializeInterior(n,z);
    
    transfer(n2,z,z1);
    
    z1[0]=z[0]=0;
    z1[n2+1]=z[n2+1]=1;
    
    
    for (i=1;i<=n2;i++)
    {
        for (j=0;j<3;j++)
        {
            Hessian[i][j]=0;
            Gtranspose_D_G[i][j]=0;
        }
        delta_z[i]=0;
    }
    
    for (i=1;i<=m;i++)
        labda[i]=w[i]=0.5;
    
    
    transfer(m,labda,labda1);
    transfer(m,w,w1);
    
    Compute_g(n,n2,z,g);
    transfer(m,g,g1);
    Compute_D(m,labda,w,D);
    
    Compute_Gtranspose_Labda(n,n2,labda,GtLabda);
    gradientInterior(n,freq,z,GtLabda,grad);
    
    Compute_Gtranspose_D_G(n,n2,D,Gtranspose_D_G);
    Compute_Hessian(n,freq,z,Gtranspose_D_G,Hessian);
    
    Compute_Gtranspose_Winverse_E(n,n2,w,Gtranspose_Winverse_E);
    Compute_Gtranspose_D_gz(n,n2,D,g,Gtranspose_D_gz);
    
    Compute_A(n2,mu,sigma,z,w,Gtranspose_D_gz,grad,Gtranspose_Winverse_E,A);
    
    // Computes expression on the right of first equation in (7.21), called A
    
    Compute_delta_z(n,n2,Hessian,A,delta_z);
    
    Compute_G_delta_z(n,n2,delta_z,G_delta_z);
    
    Compute_delta_labda(m,mu,sigma,G_delta_z,g,D,w,delta_labda);
    Compute_delta_w(m,mu,sigma,labda,w,delta_labda,D,delta_w);
    
    transfer(n2,z,z1);
    transfer(m,labda,labda1);
    transfer(m,w,w1);
    
    a=ComputeFeasibleStep(n,n2,m,z1,freq,delta_z,labda1,delta_labda,w1,delta_w);
    
    mu = DoNeighborhoodTest(m,n,n2,freq,mu,beta,z);
    
    
    j=0;
    
    while (((mu>tol) || (norm>tol) || (sum>tol)) && (j<100))
    {
        j++;
        
        mu=LoopInterior(m,n,n2,mu,sigma,beta,freq,z);
        
        phinew = phiInterior(n,n2,freq,z);
        
        sum1=0;
        for (i=1;i<=n2;i++)
            sum1=sum1+SQR(grad[i]);
        sum1=sqrt(sum1);
        
        sum2=0;
        for (i=1;i<=m;i++)
            sum2=sum2+SQR(g[i]+w[i]);
        sum2=sqrt(sum2);
        
        norm=fmax(sum1,sum2);
        
        Compute_gradient(n,freq,z,grad);
        
        sum=0;
        for (i=1;i<=n2;i++)
            sum=sum+z[i]*grad[i];
        
        gradientInterior(n,freq,z,GtLabda,grad);
        
        //printf("%10d %25.15f %25.15f %25.15f\n",j,mu,norm,phinew);
    }
    
    for (i=1;i<=n;i++)
    {
        z1[n+i]=z[2*i-1];
        z1[i]=z[2*i];
    }
    
    delete[] z1, delete[] g, delete[] g1, delete[] Gtranspose_D_gz, delete[] Gtranspose_Winverse_E;
    delete[] G_delta_z;
    
    delete[] grad, delete[] A, delete[] D, delete[] labda, delete[] w1, delete[] delta_z;
    delete[] delta_labda, delete[] delta_w, delete[] GtLabda;
    
    for (i=0;i<n2+1;i++)
        delete[] Hessian[i];
    delete[] Hessian;
    
    for (i=0;i<n2+1;i++)
        delete[] Gtranspose_D_G[i];
    delete[] Gtranspose_D_G;
}

void	Compute_delta_z(int n, int n2, double **Hessian, double A[], double delta_z[])
{
    Cholesky_sol(n2,Hessian,A,delta_z,DD,L);
}


void Compute_gradient(int n,  int **freq, double x[], double grad[])
{
    int i;
    double s;
    
    s=1.0/n;
    
    for(i=1;i<=2*n; i++)
        grad[i]=0;
    
    for(i=1;i<=n; i++)
    {
        if (freq[i][2]>0)
        {
            grad[2*i] += -s*freq[i][2]/(x[2*i]-x[2*i-1]);
            grad[2*i-1] += s*freq[i][2]/(x[2*i]-x[2*i-1]);
        }
        if (freq[i][1]>0)
            grad[2*i] += s*freq[i][1]/(1-x[2*i]);
    }
    
    for(i=1;i<=n;i++)
    {
        if (freq[i][3]>0)
        {
            if (i>1)
            {
                grad[2*i-1] += -s*freq[i][3]/(x[2*i-1]-x[2*i-3]);
                grad[2*i-3] += s*freq[i][3]/(x[2*i-1]-x[2*i-3]);
            }
            else
                grad[2*i-1] += -s*freq[i][3]/x[2*i-1];
        }
        
    }
}



void initializeInterior(int n, double z[])
{
    int i;
    
    for (i=1;i<=n;i++)
    {
        z[2*i]=i*1.0/(n+1);
        z[2*i-1]=0.9*i*1.0/(n+5);
    }
}

void Compute_Gtranspose_Labda(int n, int n2, double labda[], double GtLabda[])
{
    int i;
    
    for (i=1;i<n;i++)
        GtLabda[2*i-1] = -labda[i] + labda[i+1] + labda[n2+i];
    
    GtLabda[2*n-1] = -labda[n] + labda[n2+n];
    
    GtLabda[2] = labda[n+1] - labda[n2+1];
    
    for (i=2;i<=n;i++)
        GtLabda[2*i] = 	labda[n+i] - labda[n+i-1] - labda[n2+i];
}


void Compute_Gtranspose_D_gz(int n, int n2, double d[], double g[], double Gtranspose_D_gz[])
{
    int i;
    
    
    for (i=1;i<n;i++)
        Gtranspose_D_gz[2*i-1] = -d[i]*g[i] + d[i+1]*g[i+1] + d[n2+i]*g[n2+i];
    
    Gtranspose_D_gz[2*n-1] = -d[n]*g[n] + d[n2+n]*g[n2+n];
    
    Gtranspose_D_gz[2] = d[n+1]*g[n+1] - d[n2+1]*g[n2+1];
    
    for (i=2;i<=n;i++)
        Gtranspose_D_gz[2*i] = d[n+i]*g[n+i] - d[n+i-1]*g[n+i-1] - d[n2+i]*g[n2+i];
}

void Compute_Gtranspose_Winverse_E(int n, int n2, double w[], double Gtranspose_Winverse_E[])
{
    int i;
    
    for (i=1;i<n;i++)
        Gtranspose_Winverse_E[2*i-1] = -1/w[i] + 1/w[i+1] + 1/w[n2+i];
    
    Gtranspose_Winverse_E[2*n-1] = -1/w[n] + 1/w[n2+n];
    
    Gtranspose_Winverse_E[2] = 1/w[n+1] - 1/w[n2+1];
    
    for (i=2;i<=n;i++)
        Gtranspose_Winverse_E[2*i] = 1/w[n+i] - 1/w[n+i-1] - 1/w[n2+i];
}


void Compute_g(int n, int n2, double z[], double g[])
{
    int i;
    
    g[1] = -z[1];
    
    for (i=2;i<=n;i++)
        g[i] = z[2*i-3]-z[2*i-1];
    
    for (i=1;i<n;i++)
        g[n+i] = z[2*i]-z[2*i+2];
    
    g[n2] = z[n2]-1;
    
    for (i=1;i<=n;i++)
        g[n2+i]= z[2*i-1]-z[2*i];
}


void	Compute_G_delta_z(int n, int n2, double delta_z[], double G_delta_z[])
{
    int i;
    
    G_delta_z[1] = -delta_z[1];
    
    for (i=2;i<=n;i++)
        G_delta_z[i] = delta_z[2*i-3]-delta_z[2*i-1];
    
    for (i=1;i<n;i++)
        G_delta_z[n+i] = delta_z[2*i]-delta_z[2*i+2];
    
    G_delta_z[n2] = delta_z[n2];
    
    for (i=1;i<=n;i++)
        G_delta_z[n2+i] = delta_z[2*i-1]-delta_z[2*i];
}




void Compute_A(int n2, double mu, double sigma, double z[], double w[], double Gtranspose_D_gz[],
               double grad[], double Gtranspose_Winverse_E[], double A[])
{
    int i;
    double s;
    
    s=sigma*mu;
    
    for (i=1;i<=n2;i++)
        A[i]= -grad[i] - Gtranspose_D_gz[i] - s*Gtranspose_Winverse_E[i];
}

void Compute_delta_labda(int m, double mu, double sigma, double G_delta_z[], double g[],
                         double d[], double w[], double delta_labda[])
{
    int i;
    double s;
    
    s=sigma*mu;
    
    for (i=1;i<=m;i++)
        delta_labda[i]= d[i]*(G_delta_z[i] + g[i]) + s/w[i];
    
}

void Compute_delta_w(int m, double mu, double sigma, double labda[], double w[], double delta_labda[],
                     double D[], double delta_w[])
{
    int i;
    
    for (i=1;i<=m;i++)
        delta_w[i]=sigma*mu/labda[i]-w[i]-delta_labda[i]/D[i];
    
}



double 	phiInterior(int n, int n2, int **freq, double z[])
{
    double 	sum;
    double 	s,min;
    int 	i;
    int 	test=1;
    
    sum=0;
    min=1;
    
    i=1;
    while ((test) && (i<=n))
    {
        if (freq[i][1]>0)
        {
            s=1-z[2*i];
            if (s<min)
                min=s;
        }
        if (freq[i][2]>0)
        {
            s= z[2*i]-z[2*i-1];
            if (s<min)
                min=s;
        }
        
        if (freq[i][3]>0)
        {
            if (i>1)
                s= z[2*i-1]-z[2*i-3];
            else
                s = z[2*i-1];
            if (s<min)
                min=s;
        }
        
        if (min<0)
        {
            test=0;
            printf("constraint violation at index %10d",i);
            break;
        }
        
        i++;
    }
    
    for (i=1;i<=n; i++)
    {
        if (freq[i][1]>0)
            sum -= freq[i][1]*log(1-z[2*i]);
        if (freq[i][2]>0)
            sum -= freq[i][2]*log(z[2*i]-z[2*i-1]);
        if  (freq[i][3]>0)
        {
            if (i>1)
                sum -= freq[i][3]*log(z[2*i-1]-z[2*i-3]);
            else
                sum -= freq[i][3]*log(z[2*i-1]);
        }
    }
    
    return sum;
}


void gradientInterior(int n, int **freq, double z[], double GtLabda[], double grad[])
{
    int i;
    double s;
    
    for (i=1;i<=2*n;i++)
        grad[i]=0;
    
    s=1.0/n;
    
    for(i=1;i<=n; i++)
    {
        if (freq[i][2]>0)
        {
            grad[2*i] -= s*freq[i][2]/(z[2*i]-z[2*i-1]);
            grad[2*i-1] += s*freq[i][2]/(z[2*i]-z[2*i-1]);
        }
        if (freq[i][1]>0)
            grad[2*i] += s*freq[i][1]/(1-z[2*i]);
    }
    
    for(i=1;i<=n; i++)
    {
        if	(freq[i][3]>0)
        {
            if (i>1)
            {
                grad[2*i-1] -= s*freq[i][3]/(z[2*i-1]-z[2*i-3]);
                grad[2*i-3] += s*freq[i][3]/(z[2*i-1]-z[2*i-3]);
            }
            else
                grad[2*i-1] += -s*freq[i][3]/z[2*i-1];
        }
    }
    
    for(i=1;i<=2*n; i++)
        grad[i] = grad[i] + GtLabda[i];
}

void Compute_Hessian(int n, int **freq, double z[], double **Gtranspose_D_G, double **Hessian)
{
    int i;
    double s;
    
    s=1.0/n;
    
    for (i=1;i<=2*n;i++)
    {
        Hessian[i][0] = Gtranspose_D_G[i][0];
        Hessian[i][1] = Gtranspose_D_G[i][1];
        Hessian[i][2] = Gtranspose_D_G[i][2];
    }
    
    for(i=1;i<=n; i++)
    {
        if (freq[i][2]>0)
        {
            Hessian[2*i][0] += s*freq[i][2]/SQR(z[2*i]-z[2*i-1]);
            Hessian[2*i][1]  -=  s*freq[i][2]/SQR(z[2*i]-z[2*i-1]);
        }
        
        if (freq[i][1]>0)
            Hessian[2*i][0] += s*freq[i][1]/SQR(1-z[2*i]);
    }
    
    for(i=1;i<=n; i++)
    {
        if (freq[i][2]>0)
            Hessian[2*i-1][0]  +=  s*freq[i][2]/SQR(z[2*i]-z[2*i-1]);
        
        if (freq[i][3]>0)
        {
            if (i>1)
            {
                Hessian[2*i-1][0] += s*freq[i][3]/SQR(z[2*i-1]-z[2*i-3]);
                Hessian[2*i-3][0] += s*freq[i][3]/SQR(z[2*i-1]-z[2*i-3]);
                Hessian[2*i-1][2] -= s*freq[i][3]/SQR(z[2*i-1]-z[2*i-3]);
            }
            else
                Hessian[2*i-1][0] += s*freq[i][3]/SQR(z[2*i-1]);
        }
    }
}


void Compute_D(int m, double labda[], double w[], double d[])
{
    int i;
    for (i=1;i<=m;i++)
        d[i] = labda[i]/w[i];
}


void Compute_Gtranspose_D_G(int n, int n2, double d[], double **Gtranspose_D_G)
{
    int i;
    
    for (i=1;i<=n2;i++)
        Gtranspose_D_G[i][2]=Gtranspose_D_G[i][1]=Gtranspose_D_G[i][0]=0;
    
    for (i=1;i<n;i++)
        Gtranspose_D_G[2*i-1][0] = d[i] + d[i+1] + d[n2+i];
    
    Gtranspose_D_G[2*n-1][0] = d[n] + d[n2+n];
    
    Gtranspose_D_G[2][0] = d[n+1] + d[n2+1];
    
    for (i=2;i<=n;i++)
        Gtranspose_D_G[2*i][0] = d[n+i-1] + d[n+i] + d[n2+i];
    
    for (i=2;i<=n;i++)
    {
        Gtranspose_D_G[2*i-1][2] = -d[i];
        Gtranspose_D_G[2*i][2] = -d[n+i-1];
    }
    
    for (i=1;i<=n;i++)
        Gtranspose_D_G[2*i][1] = -d[n2+i];
    
}



double dotest(int n, int m, double z[], int **freq, double labda[], double w[])
{
    int i;
    double s,constraint;
    
    
    constraint=1.0;
    
    for (i=1;i<=n;i++)
    {
        if (freq[i][1]>0)
        {
            s=1-z[2*i];
            if (s<constraint)
                constraint = s;
        }
        if (freq[i][2]>0)
        {
            s= z[2*i]-z[2*i-1];
            if (s<constraint)
                constraint = s;
        }
        
        if (freq[i][3]>0)
        {
            if (i>1)
                s = z[2*i-1]-z[2*i-3];
            else
                s = z[2*i-1];
            if (s<constraint)
                constraint = s;
        }
        
    }
    
    
    for (i=1;i<=m;i++)
        if (fmin(labda[i],w[i])<=constraint)
            constraint=fmin(labda[i],w[i]);
    
    return constraint;
}


double 	ComputeFeasibleStep(int n, int n2, int m, double z1[], int **freq,
                            double delta_z[], double labda1[], double delta_labda[], double w1[], double delta_w[])
{
    int	i;
    double  testvalue,tol=1.0e-16;
    int	test=0;
    double  a=1.0,*z2,*labda2,*w2;
    
    z2 = new double[n2+1];
    labda2 = new double[m+1];
    w2 = new double[m+1];
    
    for (i=1;i<=n2;i++)
        z2[i]=z1[i]+delta_z[i];
    
    for (i=1;i<=m;i++)
    {
        labda2[i]=labda1[i]+delta_labda[i];
        w2[i]=w1[i]+delta_w[i];
    }
    
    while (!test)
    {
        testvalue=dotest(n,m,z2,freq,labda2,w2);
        if (testvalue>tol) test = 1;
        else
        {
            a=a/2;
            for (i=1;i<=n2;i++)
            {
                delta_z[i]=delta_z[i]/2;
                z2[i]=z1[i]+delta_z[i];
            }
            
            for (i=1;i<=m;i++)
            {
                delta_labda[i]=delta_labda[i]/2;
                delta_w[i]=delta_w[i]/2;
                labda2[i]=labda1[i]+delta_labda[i];
                w2[i]=w1[i]+delta_w[i];
            }	
        }
    }
    for (i=1;i<=n2;i++)
        z1[i]=z2[i];
    
    for (i=1;i<=m;i++)
    {
        labda1[i]=labda2[i];
        w1[i]=w2[i];
    }
    
    delete[] z2;
    delete[] labda2;
    delete[] w2;
    
    return a;
}


void transfer(int n, double a[], double b[])
{
    int	i;
    for (i = 1; i<= n;i++)	b[i] = a[i];
}



