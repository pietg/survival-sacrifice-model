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

typedef struct
{
    double x;
    int delta1;
    int delta2;
}
data;

double lambda1,lambda2;

void EM_alg(int ndata, int n, int n_mass, double vv[], int **freq, double F1[], double F2[], double f[], double *p, double **uu);
void sort_ties(int n, int delta1[], int delta2[], double xx[],
               double xx_new[], int **freq, int *n_new);
int CompareTime(const void *a, const void *b);
void sort_data(int n, double xx[], int delta1[], int delta2[]);
int  fenchelviol(int ndata, int n, int **freq, double F1[], double F2[], double f[], double tol, double *partsum, double *inprod);




// [[Rcpp::export]]

List EM(DataFrame input)
{
    int             n_mass,n,i,j,k,m2,m3,m,*delta1,*delta2;
    int             ndata,**freq;
    double          sum,*xcoor,*xx,*F2,*F1,*F1_new,*F2_new,*f,*p;
    double          *v,*w,**uu,mean1,mean2,*sum1,*sum2,*sum3;
    
    DataFrame DF = Rcpp::DataFrame(input);
    NumericVector xcoor0 = DF["V1"];
    IntegerVector delta01 = DF["V2"];
    IntegerVector delta02 = DF["V3"];
    
    // determine the sample size
    
    ndata = (int)xcoor0.size();
    
    delta1 = new int[ndata];
    delta2 = new int[ndata];
    
    xcoor = new double[ndata+2];
    
    for (i=0;i<ndata;i++)
    {
        xcoor[i]=(double)xcoor0[i];
        delta1[i]=(int)delta01[i];
        delta2[i]=(int)delta02[i];
    }
    
    sort_data(ndata,xcoor,delta1,delta2);
    
    freq = new int *[ndata];
    for (i=0;i<ndata;i++)
        freq[i]=new int[4];
    
    xx = new double[ndata];
    
    //  sort data with ties; the observation points xcoor[i] are replaced
    //  by the strictly different points xx[i]
    
    sort_ties(ndata,delta1,delta2,xcoor,xx,freq,&n);
    
    // the number of strictly different observations
    
    //printf("\nn = %5d\n\n",n);
    
    v = new double[n];
    w = new double[n];
    
    F1      = new double[n+1];
    F2      = new double[n+1];
    F1_new  = new double[n+1];
    F2_new  = new double[n+1];
    f       = new double[n+1];
    
    sum1   = new double[n+1];
    sum2 = new double[n+1];
    sum3 = new double[n+1];
    
    j=0;
    
    for (i=0;i<n;i++)
    {
        if (freq[i][3]>0)
        {
            w[j]=xx[i];
            j++;
        }
    }
    
    m3=j;
    
    j=0;
    
    for (i=0;i<n;i++)
    {
        if (freq[i][2]>0)
        {
            v[j]=xx[i];
            j++;
        }
    }
    
    m2=j;
    
    m=m2+m3;
    
    
    uu = new double *[m*m3+m2+1];
    for (i=0;i<m*m3+m2+1;i++)
        uu[i] = new double [2];
    
    k=0;
    
    for (i=0;i<m3;i++)
    {
        for (j=0;j<m2;j++)
        {
            if (v[j]<w[i])
            {
                uu[k][0]=w[i];
                uu[k][1]=v[j];
                k++;
            }
        }
        
        for (j=0;j<m3;j++)
        {
            if (w[j]<=w[i])
            {
                uu[k][0]=w[i];
                uu[k][1]=w[j];
                k++;
            }
        }
    }
    
    n_mass=k;
    
    if (freq[n-1][2]>0)
    {
        for (j=0;j<m2;j++)
        {
            uu[n_mass+j][0]=1000;
            uu[n_mass+j][1]=v[j];
        }
        n_mass += m2;
    }
    
    if (freq[n-1][1]>0)
    {
        uu[n_mass][0]=1000;
        uu[n_mass][1]=1000;
        
        n_mass++;
    }
    
    p = new double [n_mass];
    
    for (i=0;i<n_mass;i++)
        p[i]=1.0/n_mass;
    
    sum1   = new double[n];
    sum2 = new double[n];
    sum3 = new double[n];
    
    for (i=0;i<n;i++)
    {
        sum1[i]=sum2[i]=sum3[i]=0;
        for (j=0;j<n_mass;j++)
        {
            if (uu[j][0]>xx[i] && uu[j][1]>xx[i])
                sum1[i] += p[j];
            if (uu[j][0]>xx[i] && uu[j][1]<=xx[i])
                sum2[i] += p[j];
            if (uu[j][0]==xx[i])
                sum3[i] += p[j];
            
        }
        
        F1[i+1]=1-sum1[i];
        F2[i+1]=F1[i]-sum2[i];
        f[i+1]=sum3[i];
    }
    
    F1[0]=F2[0]=0;
    
    EM_alg(ndata,n,n_mass,xx,freq,F1,F2,f,p,uu);
    
    for (i=0;i<n;i++)
    {
        F1[i]=F1[i+1];
        F2[i]=F2[i+1];
        f[i]=f[i+1];
    }
    
    mean1 = (1-F1[0])*xx[0];
    mean2 = (1-F2[0])*xx[0];
    
    for (i=1;i<n;i++)
    {
        mean1 += (xx[i]-xx[i-1])*(1-F1[i]);
        mean2 += (xx[i]-xx[i-1])*(1-F2[i]);
    }
        
    sum=0;
    
    for (i=0;i<n;i++)
    {
        if (freq[i][1]>0)
            sum += freq[i][1]*log(1-F1[i]);
        if (freq[i][2]>0)
            sum += freq[i][2]*log(F1[i]-F2[i]);
        if (freq[i][3]>0 && i>0)
            sum += freq[i][3]*log(F2[i]-F2[i-1]);
    }
    
    double out2 = sum;
    
    NumericMatrix out0 = NumericMatrix(n,2);
    
    for (i=0;i<n;i++)
    {
        out0(i,0)=xx[i];
        out0(i,1)=F1[i];
    }
    
    NumericMatrix out1 = NumericMatrix(n,2);
    
    for (i=0;i<n;i++)
    {
        out1(i,0)=xx[i];
        out1(i,1)=F2[i];
    }
    
    NumericVector out3 = NumericVector(2);
    
    out3(0)=mean1;
    out3(1)=mean2;
    
    // make the list for the output, containing the two estimates and the log likelihood
    
    List EM_out = List::create(Rcpp::Named("MLE1")=out0,Rcpp::Named("MLE2")=out1,Rcpp::Named("loglikelihood")=out2,Rcpp::Named("mean")=out3);

    // free memory
    
    delete[] sum1; delete[] sum2; delete[] sum3;

    delete[] delta1; delete[] delta2; delete[] xcoor; delete[] xx;
    
    delete[] F1; delete[] F2; delete[] F1_new; delete[] F2_new; delete[] f; delete[] v; delete[] w;
    
    for (i=0;i<ndata;i++)
        delete[] freq[i];
    delete[] freq;
    
    return EM_out;
    
}

void EM_alg(int ndata, int n, int n_mass, double vv[], int **freq, double F1[], double F2[],
            double f[], double *p, double **uu)
{
    int i,j,iteration,n_It=10000;
    double sum1,sum2,sum3,sum,partialsum,inprod,tol=1.0e-10;
    
    iteration=0;
    partialsum=inprod=1;
    lambda1=lambda2=0;
    
    for (i=0;i<n;i++)
    {
      F1[i]=(i+1)*1.0/(n+1);
      F2[i]=0.9*(i+1)*1.0/(n+5);
      f[i]=0.9*1.0/(n+5);
    }
    
    //printf("iteration   partial sum     inprod\n\n");
    
    while (iteration<=n_It && fenchelviol(ndata,n,freq,F1,F2,f,tol,&partialsum,&inprod))
    {
      iteration++;
      
      for (i=0;i<n;i++)
      {
        sum1=sum2=sum3=0;
        for (j=0;j<n_mass;j++)
        {
          if (uu[j][0]>vv[i] && uu[j][1]>vv[i])
            sum1 += p[j];
          if (uu[j][0]>vv[i] && uu[j][1]<=vv[i])
            sum2 += p[j];
          if (uu[j][0]==vv[i])
            sum3 += p[j];
          
        }
        
        F1[i]=1-sum1;
        F2[i]=F1[i]-sum2;
        f[i]=sum3;
      }
      
      for (j=0;j<n_mass;j++)
      {
        sum=0;
        for (i=0;i<n;i++)
        {
          if (uu[j][0]>vv[i] && uu[j][1]>vv[i] && F1[i]<1 && freq[i][1]>0)
            sum += freq[i][1]/(1-F1[i]);
          if (uu[j][0]>vv[i] && uu[j][1]<=vv[i] && F1[i]>F2[i] && freq[i][2]>0)
            sum += freq[i][2]/(F1[i]-F2[i]);
          if (uu[j][0]==vv[i] && f[i]>0 && freq[i][3]>0)
            sum += freq[i][3]/f[i];
        }
        
        p[j]=p[j]*sum/ndata;
      }
      
      lambda1=lambda2=0;
      
      for (i=n-1;i>=0;i--)
      {
        if (freq[i][2]>0 && fabs(F1[i]-1)<1.0e-8)
          lambda1 += freq[i][2]/(ndata*(1-F2[i]));
        if (freq[i][3]>0 && fabs(F2[i]-1)<1.0e-8)
          lambda2 += freq[i][3]/(ndata*(1-F2[i-1]));
      }
    }
    
    //printf("iteration= %5d\n",iteration);
}

int  fenchelviol(int ndata, int n, int **freq, double F1[], double F2[],
                 double f[], double tol, double *partsum, double *inprod)
{
    double    sum,sum1,sum2,sum3;
    int    i,k;
    int    fenchelvioltemp=0;
    
    sum1=sum2=sum3=0;
    
    sum=lambda1;
    
    for (i=n-1;i>=0;i--)
    {
        if (freq[i][1]>0)
            sum += freq[i][1]/(ndata*(1-F1[i]));
        if (freq[i][2]>0)
            sum -= freq[i][2]/(ndata*(F1[i]-F2[i]));
        if (sum<sum1)
            sum1=sum;
    }
    
    sum = lambda1+lambda2;
    
    for (i=n-1;i>=0;i--)
    {
        if (freq[i][1]>0)
            sum += freq[i][1]/(ndata*(1-F1[i]));
        if (freq[i][2]>0)
            sum -= freq[i][2]/(ndata*(F1[i]-F2[i]));
        
        sum3=0;
        for (k=n-1;k>=i;k--)
        {
            if (freq[k][2]>0)
                sum3 += freq[k][2]/(ndata*(F1[k]-F2[k]));
            
            if (freq[k][3]>0)
                sum2 = sum3-freq[k][3]/(ndata*f[k]);
            else
                sum2=sum3;
            
            if (sum+sum2<sum1)
                sum1=sum+sum2;
        }
    }
    
    sum = lambda1+lambda2;
    
    for (i=0;i<n;i++)
    {
        if (freq[i][1]>0)
            sum += freq[i][1]*F1[i]/(ndata*(1-F1[i]));
        if (freq[i][2]>0)
            sum -= freq[i][2]*F1[i]/(ndata*(F1[i]-F2[i]));
        
        if (freq[i][2]>0)
            sum += freq[i][2]*F2[i]/(ndata*(F1[i]-F2[i]));
        
        if (freq[i][3]>0)
            sum -= freq[i][3]*F2[i]/(ndata*f[i]);
        
        if (i<n-1 && freq[i+1][3]>0)
            sum += freq[i+1][3]*F2[i]/(ndata*f[i+1]);
    }
    
    sum=fabs(sum);
    
    *inprod = sum;
    *partsum = sum1;
    
    if (sum > tol || sum1 < -tol ) fenchelvioltemp = 1;
    
    return fenchelvioltemp;
}



void sort_ties(int n, int delta1[], int delta2[], double xx[],
               double xx_new[], int **freq, int *n_new)
{
    int i,j;
    
    for (i=0;i<n;i++)
        for (j=1;j<=3;j++)
            freq[i][j]=0;
    
    if (delta1[0]==0 && delta2[0]==0)
        freq[0][1]++;
    else
    {
        if (delta1[0]==1 && delta2[0]==0)
            freq[0][2]++;
        else
            freq[0][3]++;
    }
    
    xx_new[0]=xx[0];
    
    j=0;
    
    for (i=1;i<n;i++)
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
    *n_new=j+1;
}


int CompareTime(const void *a, const void *b)
{
    if ((*(data *) a).x < (*(data *) b).x)
        return -1;
    if ((*(data *) a).x > (*(data *) b).x)
        return 1;
    return 0;
}

void sort_data(int n, double xx[], int delta1[], int delta2[])
{
    int i;
    data *obs;
    
    obs= new data[n];
    
    for (i=0;i<n;i++)
    {
        obs[i].x=xx[i];
        obs[i].delta1=delta1[i];
        obs[i].delta2=delta2[i];
    }
    
    qsort(obs,n,sizeof(data),CompareTime);
    
    for (i=0;i<n;i++)
    {
        xx[i]=obs[i].x;
        delta1[i]=obs[i].delta1;
        delta2[i]=obs[i].delta2;
    }
    
    delete[] obs;
}



