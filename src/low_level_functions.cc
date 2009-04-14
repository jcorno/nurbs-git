/* Copyright (C) 2009 Carlo de Falco
   adapted from the m-file implementation which is
   Copyright (C) 2003 Mark Spink, 2007 Daniel Claxton

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/


#include <iostream>
#include <octave/oct.h>

static int findspan(int n, int p, double u, 
		    const RowVector& U);

static void basisfun(int i, double u, int p, 
		     const RowVector& U, RowVector& N);

static double factln(int n);

static double gammaln(double xx);

static bool bspeval_bad_arguments(const octave_value_list& args);

static double bincoeff(int n, int k);

// Exports functions:
// bspeval, bspderiv, findspan

// PKG_ADD: autoload ("bspeval", "low_level_functions.oct");
DEFUN_DLD(bspeval, args, nargout,"\n\
 BSPEVAL  Evaluate B-Spline at parametric points\n\
 Calling Sequence:\n\
\n\
    p = bspeval(d,c,k,u)\n\
\n\
     INPUT:\b\
\n\
        d - Degree of the B-Spline.\n\
        c - Control Points, matrix of size (dim,nc).\n\
        k - Knot sequence, row vector of size nk.\n\
        u - Parametric evaluation points, row vector of size nu.\n\
\n\
     OUTPUT:\n\
\n\
        p - Evaluated points, matrix of size (dim,nu)\n")
{
  
  

  int       d = args(0).int_value();
  Matrix    c = args(1).matrix_value();
  RowVector k = args(2).row_vector_value();
  NDArray   u = args(3).array_value();
  
  octave_idx_type nu = u.length();
  octave_idx_type mc = c.rows(),
    nc = c.cols();

  Matrix p(mc, nu, 0.0);
  RowVector N(d+1,0.0);

  octave_value_list retval;
  if (!error_state)
    {
      if (nc + d == k.length() - 1) 
	{	 
	  int s, tmp1;
	  double tmp2;
	  
	  for (octave_idx_type col(0); col<nu; col++)
	    {	
	      s = findspan(nc-1, d, u(col), k);
	      basisfun(s, u(col), d, k, N);    
	      tmp1 = s - d;                
	      for (octave_idx_type row(0); row<mc; row++)
		{
		  double tmp2 = 0.0;
		  for ( octave_idx_type i(0); i<=d; i++)                   
		    tmp2 +=  N(i)*c(row,tmp1+i);	  
		  p(row,col) = tmp2;
		}             
	    }   
	} 
      else 
	{
	  error("inconsistent bspline data, d + columns(c) != length(k) - 1.");
	}
    }
  retval(0) = octave_value(p);
  return retval;
} 


// PKG_ADD: autoload ("bspderiv", "low_level_functions.oct");
DEFUN_DLD(bspderiv, args, nargout,"\n\
 BSPDERIV  B-Spline derivative\n     \
\n\
 Calling Sequence:\n\
\n\
          [dc,dk] = bspderiv(d,c,k)\n\
\n\
  INPUT:\n\
 \n\
    d - degree of the B-Spline\n\
    c - control points double  matrix(mc,nc)\n\
    k - knot sequence  double  vector(nk)\n\
 \n\
  OUTPUT:\n\
 \n\
    dc - control points of the derivative     double  matrix(mc,nc)\n\
    dk - knot sequence of the derivative      double  vector(nk)\n\
 \n\
  Modified version of Algorithm A3.3 from 'The NURBS BOOK' pg98.\n\
")
{
  //if (bspderiv_bad_arguments(args, nargout)) 
  //  return octave_value_list(); 
  
  int       d = args(0).int_value();
  Matrix    c = args(1).matrix_value();
  RowVector k = args(2).row_vector_value();
  octave_value_list retval;
  octave_idx_type mc = c.rows(), nc = c.cols(), nk = k.numel();
  Matrix dc (mc, nc-1, 0.0);
  RowVector dk(nk-2, 0.0);

  if (!error_state)
    {      
      double tmp;
      
      for (octave_idx_type i(0); i<=nc-2; i++)
	{
	  tmp = (double)d / (k(i+d+1) - k(i+1));
	  for ( octave_idx_type j(0); j<=mc-1; j++)
	    dc(j,i) = tmp*(c(j,i+1) - c(j,i));        
	}
      
      for ( octave_idx_type i(1); i <= nk-2; i++)
	dk(i-1) = k(i);
      
      if (nargout>1)
	retval(1) = octave_value(dk);
      retval(0) = octave_value(dc);
    }

  return(retval);
}

// Find the knot span of the parametric point u. 
//
// INPUT:
//
//   n - number of control points - 1
//   p - spline degree       
//   u - parametric point    
//   U - knot sequence
//
// RETURN:
//
//   s - knot span
//
// Algorithm A2.1 from 'The NURBS BOOK' pg68

static int findspan(int n, int p, double u, const RowVector& U)
{
  int low, high, mid;
  // special case
  if (u == U(n+1)) return(n);

  // do binary search
  low = p;
  high = n + 1;
  mid = (low + high) / 2;
  while (u < U(mid) || u >= U(mid+1))
    {

      if (u < U(mid))
	high = mid;
      else
	low = mid;
      mid = (low + high) / 2;
    }  

  return(mid);
}

// PKG_ADD: autoload ("findspan", "low_level_functions.oct");
DEFUN_DLD(findspan, args, nargout,"\n\
 FINDSPAN  Find the span of a B-Spline knot vector at a parametric point \n\
 Calling Sequence:\n\
\n\
   s = findspan(n,p,u,U)\n\
\n\
  INPUT:\n\
\n\
    n - number of control points - 1\n\
    p - spline degree\n\
    u - parametric point\n\
    U - knot sequence\n\
\n\
    U(1) <= u <= U(end)\n\
  RETURN:\n\
 \n\
    s - knot span\n\
 \n\
  Algorithm A2.1 from 'The NURBS BOOK' pg68\n")
{

  octave_value_list retval;
  int       n = args(0).idx_type_value();
  int       p = args(1).idx_type_value();
  NDArray   u = args(2).array_value();
  RowVector U = args(3).row_vector_value();
  NDArray   s(u);

  if (!error_state)
    {
      for (octave_idx_type ii(0); ii < u.length(); ii++)
	{
	  s(ii) = findspan(n, p, u(ii), U);
	}
      retval(0) = octave_value(s);
    }
  return retval;
} 


// Basis Function. 
//
// INPUT:
//
//   i - knot span  ( from FindSpan() )
//   u - parametric point
//   p - spline degree
//   U - knot sequence
//
// OUTPUT:
//
//   N - Basis functions vector[p+1]
//
// Algorithm A2.2 from 'The NURBS BOOK' pg70.

static void basisfun(int i, double u, int p, const RowVector& U, RowVector& N)
{
  int j,r;
  double saved, temp;

  // work space
  OCTAVE_LOCAL_BUFFER(double, left,  p+1);
  OCTAVE_LOCAL_BUFFER(double, right, p+1);
  
  N(0) = 1.0;
  for (j = 1; j <= p; j++)
    {
      left[j]  = u - U(i+1-j);
      right[j] = U(i+j) - u;
      saved = 0.0;
      
      for (r = 0; r < j; r++)
	{
	  temp = N(r) / (right[r+1] + left[j-r]);
	  N(r) = saved + right[r+1] * temp;
	  saved = left[j-r] * temp;
	} 
      
      N(j) = saved;
    }

}

// PKG_ADD: autoload ("basisfun", "low_level_functions.oct");
DEFUN_DLD(basisfun, args, nargout, "\n\
 Basis Function. \n\
\n\
 INPUT:\n\
\n\
   i - knot span  ( from FindSpan() )\n\
   u - parametric point\n\
   p - spline degree\n\
   U - knot sequence\n\
\n\
 OUTPUT:\n\
\n\
   N - Basis functions vector[p+1]\n\
\n\
 Algorithm A2.2 from 'The NURBS BOOK' pg70.\n")
{

  octave_value_list retval;
  NDArray   i = args(0).array_value();
  NDArray   u = args(1).array_value();
  int       p = args(2).idx_type_value();
  RowVector U = args(3).row_vector_value();
  RowVector N(p+1, 0.0);
  Matrix    B(u.length(), p+1, 0.0);
  
  if (!error_state)
    {
      for (octave_idx_type ii(0); ii < u.length(); ii++)
	{
	  basisfun(int(i(ii)), u(ii), p, U, N);
	  B.insert(N, ii, 0);
	}
      
      retval(0) = octave_value(B);
    }
  return retval;
} 

/*

%!shared n, U, p, u, s
%!test
%!  n = 3; 
%!  U = [0 0 0 1/2 1 1 1]; 
%!  p = 2; 
%!  u = linspace(0, 1, 10);  
%!  s = findspan(n, p, u, U); 
%!  assert (s, [2*ones(1, 5) 3*ones(1, 5)]);
%!test
%!  Bref = [1.00000   0.00000   0.00000
%!          0.60494   0.37037   0.02469
%!          0.30864   0.59259   0.09877
%!          0.11111   0.66667   0.22222
%!          0.01235   0.59259   0.39506
%!          0.39506   0.59259   0.01235
%!          0.22222   0.66667   0.11111
%!          0.09877   0.59259   0.30864
%!          0.02469   0.37037   0.60494
%!          0.00000   0.00000   1.00000];
%!  B = basisfun(s, u, p, U);
%!  assert (B, Bref, 1e-5);

 */

// Compute logarithm of the gamma function
// Algorithm from 'Numerical Recipes in C, 2nd Edition' pg214.
static double gammaln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6] = {76.18009172947146,-86.50532032291677,
                          24.01409824083091,-1.231739572450155,
                          0.12086650973866179e-2, -0.5395239384953e-5};
  int j;
  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x+0.5) * log(tmp);
  ser = 1.000000000190015;
  for (j=0; j<=5; j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

// computes ln(n!)
// Numerical Recipes in C
// Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215.
static double factln(int n)
{
  static int ntop = 0;
  static double a[101];
  
  if (n <= 1) return 0.0;
  while (n > ntop)
    {
      ++ntop;
      a[ntop] = gammaln(ntop+1.0);
    }
  return a[n];
}

// Computes the binomial coefficient.
//
//     ( n )      n!
//     (   ) = --------   
//     ( k )   k!(n-k)!
//
// Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215.
static double bincoeff(int n, int k)
{
  return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}


static bool bspeval_bad_arguments(const octave_value_list& args) 
{ 
  if (args.length() != 4)
    {
      error("wrong number of input arguments.");
      return true;
    }
  if (!args(0).is_real_scalar()) 
    { 
      error("degree should be a scalar."); 
      return true; 
    } 
  if (!args(1).is_real_matrix()) 
    { 
      error("the control net should be a matrix of doubles."); 
      return true; 
    } 
  if (!args(2).is_real_matrix()) 
    { 
      error("the knot vector should be a real vector."); 
      return true; 
    } 
  if (!args(3).is_real_type()) 
    { 
      error("the set of parametric points should be an array of doubles."); 
      return true; 
    } 
  return false; 
} 


