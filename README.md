# 1d-Spline
Few 1d Spline routines

The only one that is really not available elsewhere if the spline_moments() one. They are all built from the same undercomplete matrix to which is added the specific two more equations to form the specific type of splines. Note that because of taht, the natural spline function is not the most efficient possible as in principle it can be written more simply by eliminating one row and one column. There exist an iterative implementation of the spline with fix first derivatives, should be added in the future and compared to the matrix implementation.
