#include <stdio.h>
#include <stdlib.h>

void spline_coefficients_deriv(const double *xx,const double *yy,const double *fp,double *vec_coeff,const int size)
{

	//Spline with given first derivatives
	double FP0=fp[0];
	double FPN=fp[1];

	int n=size-1;

	double *h,*alpha,*mu,*z,*l;

	h=(double *) malloc(n*sizeof(double));
	alpha=(double *) malloc((n+1)*sizeof(double));
	mu=(double *) malloc(n*sizeof(double));
	z=(double *) malloc((n+1)*sizeof(double));
	l=(double *) malloc((n+1)*sizeof(double));

	int r;
	for (r=0;r<n;r++)
	{
		h[r]=xx[r+1]-xx[r];
		if (r==0)
  			alpha[0]=3.0*(yy[1]-yy[0])/h[0] - 3.0*FP0;
		else
  			alpha[r]=3.0/h[r]*(yy[r+1]-yy[r]) - 3.0/h[r-1]*(yy[r]-yy[r-1] );
	}
	alpha[n]=3.0*FPN - 3.0*( yy[n] - yy[n-1] )/h[n-1];

	l[0]=2*h[0];
	mu[0]=0.5;
	z[0]=alpha[0]/l[0];

	for (r=1;r<n;r++)
	{
    		l[r]=2.0*(xx[r+1]-xx[r-1])-h[r-1]*mu[r-1];
    		mu[r]=h[r]/l[r];
    		z[r]=(alpha[r]-h[r-1]*z[r-1]  )/l[r];
	}
	l[n]=h[n-1]*(2.0-mu[n-1]);
	z[n]=( alpha[n] - h[n-1]*z[n-1]  )/l[n];

	double *c,*b,*d;
	c=(double *) malloc((n+1)*sizeof(double));
	b=(double *) malloc(n*sizeof(double));
	d=(double *) malloc(n*sizeof(double));

	c[n]=z[n];
    	for (r=0;r<n;r++)
	{
    		c[(n-1)-r]=z[(n-1)-r] - mu[(n-1)-r]*c[n-r];
    		b[(n-1)-r]=(yy[n-r]-yy[(n-1)-r])/h[(n-1)-r] - h[(n-1)-r]*( c[n-r]+2.0*c[(n-1)-r] )/3.0;
    		d[(n-1)-r]=( c[n-r] - c[(n-1)-r] )/(3.0*h[(n-1)-r]);
	}

	for (r=0;r<n;r++)
	{
   		vec_coeff[4*r]=yy[r];
   		vec_coeff[4*r+1]=b[r];
   		vec_coeff[4*r+2]=c[r];
   		vec_coeff[4*r+3]=d[r];
	}

	free(h);
	free(alpha);
	free(mu);
	free(z);
	free(l);
	free(c);
	free(b);
	free(d);

}
