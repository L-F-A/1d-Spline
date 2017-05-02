# 1d-Spline

Few 1d Spline routines

The only one that is really not available elsewhere if the spline_moments() one. 

      1. The methods are all built from the same undercomplete matrix to which is added the specific two more equations
         to form the specific type of splines. Note that because of that, the natural spline function is not the most 
         efficient possible as in principle it can be written more simply by eliminating one row and one column. 

      2. For the case with fixed first derivatives, there exist an iterative algorithm. It is possible to use it by 
         choosing the parameter iterative=True in the function spline_1stderiv.
            
            -The core function for the iterative algorithm is written in C. To use it, the C codes must be compiled in the same directory                where spline_lib.py and firstDerivSpline.c are. The compilation line is just for Linux (I guess Mac to):

                                    gcc -shared -fPIC firstDerivSpline.c -o firstDerivSpline.so

Then evrerything should be ok.

-The function locate() from Standard Useful repository needs to be added as it is loaded by the file.
