import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg.dsolve import linsolve
import math
from StandardUseful import locate

####################################################################################################
#         This is a simple library that offers few cublic splines for 1d interpolation		   #
#												   #
#           Louis-Francois Arsenault, Columbia University la2518@columbia.edu (2013-2017)	   #
####################################################################################################
#												   #
#	1. Natural as spline_natural(x,f,n)							   #
#												   #
#	2. End points derivatives fixed as spline_1stderiv(x,f,df0,df1,n)			   #
#												   #
#	3. Parabolically-terminated as spline_parabolic(x,f,n)					   #
#												   #
#	4. Special form where one knows two constants such that f'(x0)+f'(xf) = mu1 and		   #
#	   -f''(x0)-f''(xf) = mu2 as spline_moments(x,f,mu1,mu2,n), see for example Appendix B of  #
# 	   http://e-collection.library.ethz.ch/eserv/eth:31103/eth-31103-02.pdf			   #
#												   #
#	INPUTS:											   #
#	  											   #
#	  x 	  : the points where the function is known					   #
#												   #
#	  f 	  : the values of the function a the x						   #
#												   #
#	  n 	  : how many points								   #
#	Possible:										   #
#												   #
#	  df0 	  : first derivative at x0							   #
#												   #
#	  df1 	  : first derivate at xf							   #
#												   #
#	  mu1,mu2 : knowns constants								   #
####################################################################################################

def UnderCompletedAmatrix(n):
    A = lil_matrix((n,n))
    Ld = 4.*np.ones(n)
    L = np.ones(n-1)
    A.setdiag(Ld)
    A.setdiag(L, k=1)
    A[0,0] = 0.
    A[n-1,n-1] = 0.    
    A[0,1] = 0.
    A.setdiag(L, k=-1)    
    return A

def UndercompletedBmatrix(n,f):
    B = np.zeros(n)
    for r in range(1,n-1):
        B[r] = f[r+1] - 2*f[r] + f[r-1]
    return B

def spline_moments(x,f,mu1,mu2,n):
#This is a spline for which one knows two constants mu1 and mu2 such that
#the function has properties f'(x0)+f'(xf) = mu1 and -f''(x0)-f''(xf) = mu2
#and thus we can close the system of equations of the spline using mu1 and mu2
#
#Example: This is useful in quantum many-body physics where f(x) is G(tau) the Green's 
#function in imaginary time and mu1 and mu2 are the first and second moments 
#(counting starts at 0) 
 
    A = UnderCompletedAmatrix(n)
    B = UndercompletedBmatrix(n,f)     
    
    h = (x[-1]-x[0])/(n-1);

    B[0] = -mu2*h
    B[-1] = -( f[1]-f[0] + f[-1] - f[-2]  ) + mu1*h
    B = 6./h/h*B

    A[0,0] = 6/h
    A[0,n-1] = 6/h
    A[n-1,0] = -2.
    A[n-1,1] = -1.
    A[n-1,n-1] = 2.
    A[n-1,n-2] = 1.

    A = A.tocsr()
    M = linsolve.spsolve(A, B)

    a = np.zeros(n-1)
    b = np.zeros(n-1)
    c = np.zeros(n-1)
    d = f[0:n-1].copy()
    
    for r in range(0,n-1):
        a[r] = (M[r+1]-M[r])/6./h;
        b[r] = 0.5*M[r];
        c[r] = ( f[r+1]-f[r] )/h - h/6.*(M[r+1]+2.*M[r])

    COEFFS = np.zeros((n-1,4))
    COEFFS[:,0] = a.copy()
    COEFFS[:,1] = b.copy() 
    COEFFS[:,2] = c.copy()
    COEFFS[:,3] = d.copy()

    return COEFFS
def spline_natural(x,f,n):
#Second derivatives set to zero at both end   

    A = UnderCompletedAmatrix(n)
    B = UndercompletedBmatrix(n,f)     
    
    h = (x[-1]-x[0])/(n-1);

    B = 6./h/h*B
    A[0,0] = 1.
    A[n-1,n-1] = 1.
    
    A = A.tocsr()
    M = linsolve.spsolve(A,B)
    
    a = np.zeros(n-1)
    b = np.zeros(n-1)
    c = np.zeros(n-1)
    d = f[0:n-1].copy()
    
    for r in range(0,n-1):
        a[r] = (M[r+1]-M[r])/6./h;
        b[r] = 0.5*M[r];
        c[r] = ( f[r+1]-f[r] )/h - h/6.*(M[r+1]+2.*M[r])

    
    COEFFS = np.zeros((n-1,4))
    COEFFS[:,0] = a.copy()
    COEFFS[:,1] = b.copy() 
    COEFFS[:,2] = c.copy()
    COEFFS[:,3] = d.copy()

    return COEFFS

def spline_parabolic(x,f,n):
#parabolically-terminated cubic spline   

    A = UnderCompletedAmatrix(n)
    B = UndercompletedBmatrix(n,f)     
    
    h = (x[-1]-x[0])/(n-1);

    B = 6./h/h*B
    A[0,0] = 1.
    A[0,1] = -1.
    A[n-1,n-1] = 1.
    A[n-1,n-2] = -1.
    
    A = A.tocsr()
    M = linsolve.spsolve(A,B)
    
    a = np.zeros(n-1)
    b = np.zeros(n-1)
    c = np.zeros(n-1)
    d = f[0:n-1].copy()
    
    for r in range(0,n-1):
        a[r] = (M[r+1]-M[r])/6./h;
        b[r] = 0.5*M[r];
        c[r] = ( f[r+1]-f[r] )/h - h/6.*(M[r+1]+2.*M[r])

    
    COEFFS = np.zeros((n-1,4))
    COEFFS[:,0] = a.copy()
    COEFFS[:,1] = b.copy() 
    COEFFS[:,2] = c.copy()
    COEFFS[:,3] = d.copy()

    return COEFFS    

def spline_1stderiv(x,f,df0,df1,n):
#spline where first derivatives at both end are fixed

    #df0 is the value of the first derivative at the beginning
    #df1 is the value of the first derivative at the end
    
    A = UnderCompletedAmatrix(n)
    B = UndercompletedBmatrix(n,f)     
    
    h = (x[-1]-x[0])/(n-1);

    B[0] = f[1] - f[0] - df0*h
    B[-1] = df1*h - f[-1] + f[-2] 
    B = 6./h/h*B

    A[0,0] = 2.
    A[0,1] = 1.
    A[n-1,n-1] = 2.
    A[n-1,n-2] = 1.

    A = A.tocsr()
    M = linsolve.spsolve(A, B)

    a = np.zeros(n-1)
    b = np.zeros(n-1)
    c = np.zeros(n-1)
    d = f[0:n-1].copy()
    
    for r in range(0,n-1):
        a[r] = (M[r+1]-M[r])/6./h;
        b[r] = 0.5*M[r];
        c[r] = ( f[r+1]-f[r] )/h - h/6.*(M[r+1]+2.*M[r])

    
    COEFFS = np.zeros((n-1,4))
    COEFFS[:,0] = a.copy()
    COEFFS[:,1] = b.copy() 
    COEFFS[:,2] = c.copy()
    COEFFS[:,3] = d.copy()

    return COEFFS

def spline_eval(x,COEFFS,xx):
#Evaluates the spline prediction at point x
    y = np.zeros(len(x))
    for r in range(0,len(x)):
        ind = locate(x[r],xx)
	if ind[0] == len(xx)-1:
		ind[0] = ind[0]-1
        a = COEFFS[ind[0],0]
        b = COEFFS[ind[0],1]
        c = COEFFS[ind[0],2]
        d = COEFFS[ind[0],3]
        y[r] = a*( x[r] - xx[ind[0]]  )**3 + b*(x[r]-xx[ind[0]])**2 + c*(x[r]-xx[ind[0]]) + d
    return y

