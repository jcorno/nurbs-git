/* Copyright (C) 2009 Carlo de Falco
   some functions are adapted from the m-file implementation which is
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
#include <octave/oct-map.h>
#include <octave/parse.h>

static int findspan(int n, int p, double u, const RowVector& U);
static void basisfun(int i, double u, int p, const RowVector& U, RowVector& N);
static void basisfunder (int i, int pl, double uu, const RowVector& u_knotl, int nders, NDArray& dersv);
static double factln(int n);
static double gammaln(double xx);
static double bincoeff(int n, int k);
static bool bspeval_bad_arguments(const octave_value_list& args);

// Exported functions:
// bspeval, bspderiv, findspan, basisfun, __nrb_srf_basisfun__, __nrb_srf_basisfun_der__

// PKG_ADD: autoload ("basisfunder", "low_level_functions.oct");
DEFUN_DLD(basisfunder, args, nargout,"\n\
 BASISFUNDER:  B-Spline Basis function derivatives\n\
\n\
 Calling Sequence:\n\
\n\
   ders = basisfunder (i, pl, u, k, nd)\n\
\n\
    INPUT:\n\
\n\
      i   - knot span\n\
      pl  - degree of curve\n\
      u   - parametric points\n\
      k   - knot vector\n\
      nd  - number of derivatives to compute\n\
\n\
    OUTPUT:\n\
\n\
      ders - ders(n, i, :) (i-1)-th derivative at n-th point\n\
")
{
  octave_value_list retval;

  const RowVector i = args(0).row_vector_value ();
  int pl = args(1).int_value ();
  const RowVector u = args(2).row_vector_value ();
  const RowVector U = args(3).row_vector_value ();
  int nd = args(4).int_value ();

  if (!error_state)
    {
      if (i.length () != u.length ())
	print_usage ();
 
      NDArray dersv (dim_vector (i.length (), nd+1, pl+1), 0.0);
      NDArray ders(dim_vector(nd+1, pl+1), 0.0);
      for ( octave_idx_type jj(0); jj < i.length (); jj++)
	{
	  basisfunder (int (i(jj)), pl, u(jj), U, nd, ders);

	  for (octave_idx_type kk(0); kk < nd+1; kk++)
	    for (octave_idx_type ll(0); ll < pl+1; ll++)
	      {
		dersv(jj, kk, ll) = ders(kk, ll);
	      }
	}
      retval(0) = dersv;
    }
  return retval;
}


// PKG_ADD: autoload ("__nrb_srf_basisfun_der__", "low_level_functions.oct");
DEFUN_DLD(__nrb_srf_basisfun_der__, args, nargout,"\
 __NRB_SRF_BASISFUN_DER__:  Undocumented private function	\
")
{
  //function [Bu, Bv, N] = __nrb_srf_basisfun_der__ (points, nrb);

  octave_value_list retval, newargs;

  const NDArray points = args(0).array_value();
  const Octave_map nrb = args(1).map_value();

  if (!error_state) 
    {
      const Cell knots = nrb.contents("knots")(0).cell_value();
      const NDArray coefs = nrb.contents("coefs")(0).array_value();
      octave_idx_type m   = (nrb.contents("number")(0).vector_value())(0) - 1; // m    = size (nrb.coefs, 2) -1;
      octave_idx_type n   = (nrb.contents("number")(0).vector_value())(1) - 1; // n    = size (nrb.coefs, 3) -1;
      octave_idx_type p = (nrb.contents("order")(0).vector_value())(0) - 1;    // p    = nrb.order(1) -1;
      octave_idx_type q = (nrb.contents("order")(0).vector_value())(1) - 1;    // q    = nrb.order(2) -1;

      Array<idx_vector> idx(2, idx_vector(':')); 
      idx(0) = 1;
      const NDArray u(points.index (idx).squeeze ()); // u = points(1,:);
      
      idx(0) = 2;
      const NDArray v(points.index (idx).squeeze ()); // v = points(2,:);      
      
      octave_idx_type npt = u.length (); // npt = length(u);

      RowVector M(p+1, 0.0), N (q+1, 0.0);
      Matrix Nout(npt, (p+1)*(q+1), 0.0);
      Matrix Bu(npt, (p+1)*(q+1), 0.0);
      Matrix Bv(npt, (p+1)*(q+1), 0.0);
      RowVector Denom(npt, 0.0);
      RowVector Denom_du(npt, 0.0);
      RowVector Denom_dv(npt, 0.0);
      Matrix Num(npt, (p+1)*(q+1), 0.0);
      Matrix Num_du(npt, (p+1)*(q+1), 0.0);
      Matrix Num_dv(npt, (p+1)*(q+1), 0.0);

      const RowVector U(knots(0).row_vector_value ()); // U = nrb.knots{1};

      const RowVector V(knots(1).row_vector_value ()); // V = nrb.knots{2};
      
      Array<idx_vector> idx2(3, idx_vector(':')); idx2(0) = 4;
      NDArray w (coefs.index (idx2).squeeze ()); // w = squeeze(nrb.coefs(4,:,:));
      
      RowVector spu(u);
      for (octave_idx_type ii(0); ii < npt; ii++)
	{
	  spu(ii) = findspan(m, p, u(ii), U);
	} // spu  =  findspan (m, p, u, U); 

      newargs(3) = U; newargs(2) = p; newargs(1) = u; newargs(0) = spu;
      Matrix Ik = feval (std::string("numbasisfun"), newargs, 1)(0).matrix_value (); // Ik = numbasisfun (spu, u, p, U);

      RowVector spv(v);
      for (octave_idx_type ii(0); ii < v.length(); ii++)
	{
	  spv(ii) = findspan(n, q, v(ii), V);
	} // spv  =  findspan (n, q, v, V);

      newargs(3) = V; newargs(2) = q; newargs(1) = v; newargs(0) = spv;
      Matrix Jk = feval (std::string("numbasisfun"), newargs, 1)(0).matrix_value (); // Jk = numbasisfun (spv, v, q, V);

      Matrix NuIkuk(npt, p+1, 0.0);
      for (octave_idx_type ii(0); ii < npt; ii++)
	{
	  basisfun (int(spu(ii)), u(ii), p, U, M);
	  NuIkuk.insert (M, ii, 0);
	} // NuIkuk = basisfun (spu, u, p, U);

      Matrix NvJkvk(v.length (), q+1, 0.0);
      for (octave_idx_type ii(0); ii < npt; ii++)
	{
	  basisfun(int(spv(ii)), v(ii), q, V, N);
	  NvJkvk.insert (N, ii, 0);
	} // NvJkvk = basisfun (spv, v, q, V);

     
      newargs(4) = 1; newargs(3) = U; newargs(2) = u; newargs(1) = p; newargs(0) = spu;
      NDArray NuIkukprime_tmp = feval (std::string("basisfunder"), newargs, 1)(0).array_value (); //   NuIkukprime = basisfunder (spu, p, u, U, 1);
      idx2(0) = idx_vector(':'); idx2(1) = 2; idx2(2) = idx_vector(':');
      NDArray NuIkukprime (NuIkukprime_tmp.index (idx2).squeeze ()); // NuIkukprime = squeeze(NuJkukprime(:,2,:));
      
      newargs(4) = 1; newargs(3) = V; newargs(2) = v; newargs(1) = q; newargs(0) = spv;
      NDArray NvJkvkprime_tmp = feval (std::string("basisfunder"), newargs, 1)(0).array_value (); //   NvJkvkprime = basisfunder (spv, q, v, V, 1);
      idx2(0) = idx_vector(':'); idx2(1) = 2; idx2(2) = idx_vector(':');
      NDArray NvJkvkprime (NvJkvkprime_tmp.index (idx2).squeeze ()); // NvJkvkprime = squeeze(NvJkvkprime(:,2,:));
      
      for (octave_idx_type k(0); k < npt; k++) 
	for (octave_idx_type ii(0); ii < p+1; ii++) 
	  for (octave_idx_type jj(0); jj < q+1; jj++) 
	    {
	      Num(k, ii+jj*(p+1)) = NuIkuk(k, ii) * NvJkvk(k, jj) * w(Ik(k, ii), Jk(k, jj));
	      Denom(k) += Num(k, ii+jj*(p+1));

	      Num_du(k, ii+jj*(p+1)) = NuIkukprime(k, ii) * NvJkvk(k, jj) * w(Ik(k, ii), Jk(k, jj));
	      Denom_du(k) += Num_du(k, ii+jj*(p+1));

	      Num_dv(k, ii+jj*(p+1)) = NuIkuk(k, ii) * NvJkvkprime(k, jj) * w(Ik(k, ii), Jk(k, jj));
	      Denom_dv(k) += Num_dv(k, ii+jj*(p+1));
	    }

      for (octave_idx_type k(0); k < npt; k++) 
	for (octave_idx_type ii(0); ii < p+1; ii++) 
	  for (octave_idx_type jj(0); jj < q+1; jj++) 
	    {
	      Bu(k, octave_idx_type(ii+(p+1)*jj))  = (Num_du(k, ii+jj*(p+1))/Denom(k) - Denom_du(k)*Num(k, ii+jj*(p+1))/(Denom(k)*Denom(k))); 
	      Bv(k, octave_idx_type(ii+(p+1)*jj))  = (Num_dv(k, ii+jj*(p+1))/Denom(k) - Denom_dv(k)*Num(k, ii+jj*(p+1))/(Denom(k)*Denom(k))); 
	      Nout(k, octave_idx_type(ii+(p+1)*jj))= Ik(k, ii)+(m+1)*Jk(k, jj)+1;
	    }

      //   for k=1:npt
      //     [Ika, Jkb] = meshgrid(Ik(k, :), Jk(k, :)); 
      
      //     N(k, :)    = sub2ind([m+1, n+1], Ika(:)+1, Jkb(:)+1);
      //     wIkaJkb(1:p+1, 1:q+1) = reshape (w(N(k, :)), p+1, q+1); 
      
      //     Num    = (NuIkuk(k, :).' * NvJkvk(k, :)) .* wIkaJkb;
      //     Num_du = (NuIkukprime(k, :).' * NvJkvk(k, :)) .* wIkaJkb;
      //     Num_dv = (NuIkuk(k, :).' * NvJkvkprime(k, :)) .* wIkaJkb;
      //     Denom  = sum(sum(Num));
      //     Denom_du = sum(sum(Num_du));
      //     Denom_dv = sum(sum(Num_dv));
      
      //     Bu(k, :) = (Num_du/Denom - Denom_du.*Num/Denom.^2)(:).';
      //     Bv(k, :) = (Num_dv/Denom - Denom_dv.*Num/Denom.^2)(:).';
      //   end
      

      retval(2) = Nout;
      retval(1) = Bv;
      retval(0) = Bu;

    }
  return retval;
}

// PKG_ADD: autoload ("__nrb_srf_basisfun__", "low_level_functions.oct");
DEFUN_DLD(__nrb_srf_basisfun__, args, nargout,"\
 __NRB_SRF_BASISFUN__:  Undocumented private function\
")
{

  octave_value_list retval, newargs;

  const NDArray points = args(0).array_value();
  const Octave_map nrb = args(1).map_value();

  if (!error_state) 
    {

      const Cell knots = nrb.contents("knots")(0).cell_value();
      const NDArray coefs = nrb.contents("coefs")(0).array_value();
      octave_idx_type m   = (nrb.contents("number")(0).vector_value())(0) - 1; // m    = size (nrb.coefs, 2) -1;
      octave_idx_type n   = (nrb.contents("number")(0).vector_value())(1) - 1; // n    = size (nrb.coefs, 3) -1;
      octave_idx_type p = (nrb.contents("order")(0).vector_value())(0) - 1;    // p    = nrb.order(1) -1;
      octave_idx_type q = (nrb.contents("order")(0).vector_value())(1) - 1;    // q    = nrb.order(2) -1;

      Array<idx_vector> idx(2, idx_vector(':')); 
      idx(0) = 1;
      const NDArray u(points.index (idx).squeeze ()); // u = points(1,:);

      idx(0) = 2;
      const NDArray v(points.index (idx).squeeze ()); // v = points(2,:);      

      octave_idx_type npt = u.length (); // npt = length(u);
      RowVector M(p+1, 0.0), N (q+1, 0.0);
      Matrix RIkJk(npt, (p+1)*(q+1), 0.0);
      Matrix indIkJk(npt, (p+1)*(q+1), 0.0);
      RowVector denom(npt, 0.0);

      const RowVector U(knots(0).row_vector_value ()); // U = nrb.knots{1};

      const RowVector V(knots(1).row_vector_value ()); // V = nrb.knots{2};
      
      Array<idx_vector> idx2(3, idx_vector(':')); idx2(0) = 4;
      NDArray w (coefs.index (idx2).squeeze ()); // w = squeeze(nrb.coefs(4,:,:));
      
      RowVector spu(u);
      for (octave_idx_type ii(0); ii < npt; ii++)
	{
	  spu(ii) = findspan(m, p, u(ii), U);
	} // spu  =  findspan (m, p, u, U); 

      newargs(3) = U; newargs(2) = p; newargs(1) = u; newargs(0) = spu;
      Matrix Ik = feval (std::string("numbasisfun"), newargs, 1)(0).matrix_value (); // Ik = numbasisfun (spu, u, p, U);

      RowVector spv(v);
      for (octave_idx_type ii(0); ii < v.length(); ii++)
	{
	  spv(ii) = findspan(n, q, v(ii), V);
	} // spv  =  findspan (n, q, v, V);

      newargs(3) = V; newargs(2) = q; newargs(1) = v; newargs(0) = spv;
      Matrix Jk = feval (std::string("numbasisfun"), newargs, 1)(0).matrix_value (); // Jk = numbasisfun (spv, v, q, V);

      Matrix NuIkuk(npt, p+1, 0.0);
      for (octave_idx_type ii(0); ii < npt; ii++)
	{
	  basisfun (int(spu(ii)), u(ii), p, U, M);
	  NuIkuk.insert (M, ii, 0);
	} // NuIkuk = basisfun (spu, u, p, U);

      Matrix NvJkvk(v.length (), q+1, 0.0);
      for (octave_idx_type ii(0); ii < npt; ii++)
	{
	  basisfun(int(spv(ii)), v(ii), q, V, N);
	  NvJkvk.insert (N, ii, 0);
	} // NvJkvk = basisfun (spv, v, q, V);


      for (octave_idx_type k(0); k < npt; k++) 
	for (octave_idx_type ii(0); ii < p+1; ii++) 
	  for (octave_idx_type jj(0); jj < q+1; jj++) 
	    denom(k) += NuIkuk(k, ii) * NvJkvk(k, jj) * w(Ik(k, ii), Jk(k, jj));

      
      for (octave_idx_type k(0); k < npt; k++) 
	for (octave_idx_type ii(0); ii < p+1; ii++) 
	  for (octave_idx_type jj(0); jj < q+1; jj++) 
	    {
	      RIkJk(k, octave_idx_type(ii+(p+1)*jj))  = NuIkuk(k, ii)*NvJkvk(k, jj) * w(Ik(k, ii), Jk(k, jj))/denom(k); 
	      indIkJk(k, octave_idx_type(ii+(p+1)*jj))= Ik(k, ii)+(m+1)*Jk(k, jj)+1;
	    }

      // for k=1:npt
      //       [Jkb, Ika] = meshgrid(Jk(k, :), Ik(k, :)); 
      //       indIkJk(k, :)    = sub2ind([m+1, n+1], Ika(:)+1, Jkb(:)+1);
      //       wIkaJkb(1:p+1, 1:q+1) = reshape (w(indIkJk(k, :)), p+1, q+1); 
      
      //       NuIkukaNvJkvk(1:p+1, 1:q+1) = (NuIkuk(k, :).' * NvJkvk(k, :));
      //       RIkJk(k, :) = (NuIkukaNvJkvk .* wIkaJkb ./ sum(sum(NuIkukaNvJkvk .* wIkaJkb)))(:).';
      //     end
      
      retval(0) = RIkJk; // B = RIkJk;
      retval(1) = indIkJk; // N = indIkJk;

    }
  return retval;
}


// PKG_ADD: autoload ("bspeval", "low_level_functions.oct");
DEFUN_DLD(bspeval, args, nargout,"\
 BSPEVAL:  Evaluate B-Spline at parametric points\n\
\n\
\n\
 Calling Sequence:\n					\
\n							\
   p = bspeval(d,c,k,u)\n				\
\n							\
    INPUT:\n						\
\n							\
       d - Degree of the B-Spline.\n			\
       c - Control Points, matrix of size (dim,nc).\n	\
       k - Knot sequence, row vector of size nk.\n			\
       u - Parametric evaluation points, row vector of size nu.\n	\
 \n									\
    OUTPUT:\n								\
\n									\
       p - Evaluated points, matrix of size (dim,nu)\n			\
")
{
  
  

  int       d = args(0).int_value();
  const Matrix    c = args(1).matrix_value();
  const RowVector k = args(2).row_vector_value();
  const NDArray   u = args(3).array_value();
  
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
 BSPDERIV:  B-Spline derivative\n\
\n\
\n\
 Calling Sequence:\n				\
\n						\
          [dc,dk] = bspderiv(d,c,k)\n		\
\n						\
  INPUT:\n					\
 \n						\
    d - degree of the B-Spline\n		\
    c - control points double  matrix(mc,nc)\n	\
    k - knot sequence  double  vector(nk)\n	\
 \n						\
  OUTPUT:\n					\
 \n									\
    dc - control points of the derivative     double  matrix(mc,nc)\n	\
    dk - knot sequence of the derivative      double  vector(nk)\n	\
 \n									\
  Modified version of Algorithm A3.3 from 'The NURBS BOOK' pg98.\n	\
")
{
  //if (bspderiv_bad_arguments(args, nargout)) 
  //  return octave_value_list(); 
  
  int       d = args(0).int_value();
  const Matrix    c = args(1).matrix_value();
  const RowVector k = args(2).row_vector_value();
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


static int findspan(int n, int p, double u, const RowVector& U)

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
DEFUN_DLD(findspan, args, nargout,"\
FINDSPAN: Find the span of a B-Spline knot vector at a parametric point\n \
\n\
\n\
Calling Sequence:\n							\
\n									\
   s = findspan(n,p,u,U)\n						\
\n									\
  INPUT:\n								\
\n									\
    n - number of control points - 1\n					\
    p - spline degree\n							\
    u - parametric point\n						\
    U - knot sequence\n							\
\n									\
    U(1) <= u <= U(end)\n						\
  RETURN:\n								\
 \n									\
    s - knot span\n							\
 \n									\
  Algorithm A2.1 from 'The NURBS BOOK' pg68\n				\
")
{

  octave_value_list retval;
  int       n = args(0).idx_type_value();
  int       p = args(1).idx_type_value();
  const NDArray   u = args(2).array_value();
  const RowVector U = args(3).row_vector_value();
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



static void basisfun(int i, double u, int p, const RowVector& U, RowVector& N)
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
 BASISFUN: Compute B-Spline Basis Functions \n	\
\n						\
 INPUT:\n					\
\n						\
   i - knot span  ( from FindSpan() )\n		\
   u - parametric point\n			\
   p - spline degree\n				\
   U - knot sequence\n				\
\n						\
 OUTPUT:\n					\
\n						\
   N - Basis functions vector[p+1]\n		\
\n						\
 Algorithm A2.2 from 'The NURBS BOOK' pg70.\n	\
")
{

  octave_value_list retval;
  const NDArray   i = args(0).array_value();
  const NDArray   u = args(1).array_value();
  int       p = args(2).idx_type_value();
  const RowVector U = args(3).row_vector_value();
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


static double gammaln(double xx)
// Compute logarithm of the gamma function
// Algorithm from 'Numerical Recipes in C, 2nd Edition' pg214.
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

static double factln(int n)
// computes ln(n!)
// Numerical Recipes in C
// Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215.
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

static double bincoeff(int n, int k)
// Computes the binomial coefficient.
//
//     ( n )      n!
//     (   ) = --------   
//     ( k )   k!(n-k)!
//
// Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215.
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


static void basisfunder (int i, int pl, double u, const RowVector& u_knotl, int nders, NDArray& ders)
{
  
  //     ders = zeros(nders+1,pl+1);
  Matrix ndu(octave_idx_type(pl+1), octave_idx_type(pl+1), 0.0); // ndu = zeros(pl+1,pl+1);
  RowVector left(octave_idx_type(pl+1), 0.0);                    // left = zeros(pl+1);
  RowVector right(left);                                         // right = zeros(pl+1);
  Matrix a(2, octave_idx_type(pl+1), 0.0);                       // a = zeros(2,pl+1);
  double saved = 0.0, d = 0.0, temp = 1.0;
  octave_idx_type s1(0), s2(1), rk, pk, j, k, r, j1, j2;

  
  ndu(0,0) = 1;                                                   // ndu(1,1) = 1;
  
  for (j=1; j<=pl; j++)                                           // for j = 1:pl
    {
      left(j) = u - u_knotl(i+1-j);                               // left(j+1) = u - u_knotl(i+1-j);
      right(j) = u_knotl(i+j) - u;                                // right(j+1) = u_knotl(i+j) - u;
      saved = 0.0;                                                // saved = 0;
      for (r=0; r<=j-1; r++)                                      // for r = 0:j-1
	{
	  ndu(j, r) = right(r+1) + left(j-r);                     // ndu(j+1,r+1) = right(r+2) + left(j-r+1);
	  temp = ndu(r,j-1)/ndu(j,r);                             // temp = ndu(r+1,j)/ndu(j+1,r+1);
	  ndu(r,j) = saved + right(r+1)*temp;                     // ndu(r+1,j+1) = saved + right(r+2)*temp;
	  saved = left(j-r)*temp;                                 // saved = left(j-r+1)*temp;
	}                                                         // end
      ndu(j,j) = saved;                                           // ndu(j+1,j+1) = saved;
    }                                                             // end

  for (j=0; j<=pl; j++)                                           // for j = 0:pl
    ders(0,j) = ndu(j,pl);                                        // ders(1,j+1) = ndu(j+1,pl+1);
                                                                  // end
      
  for (r=0; r<=pl; r++)                                           // for r = 0:pl
    {
      s1 = 0;                                                     // s1 = 0;
      s2 = 1;                                                     // s2 = 1;
      a(0,0) = 1;                                                 // a(1,1) = 1;
        for (k=1; k<=nders; k++)                                  // for k = 1:nders %compute kth derivative
	  {

	    d = 0.0;                                              // d = 0;
	    rk = r-k;                                             // rk = r-k;
	    pk = pl - k;                                          // pk = pl-k;
	    
	    if (r >= k)                                           // if (r >= k)
	      {
		a(s2, 0) = a(s1, 0)/ndu(pk+1,rk);                 // a(s2+1,1) = a(s1+1,1)/ndu(pk+2,rk+1);
		d = a(s2, 0)*ndu(rk,pk);                          // d = a(s2+1,1)*ndu(rk+1,pk+1);
	      }                                                   // end
	    
	    if (rk >= -1)                                         // if (rk >= -1)
	      j1 = 1;                                             // j1 = 1;
	    else                                                  // else 
	      j1 = -rk;                                           // j1 = -rk;
	                                                          // end

	    if ((r-1) <= pk)                                      // if ((r-1) <= pk)
	      j2 = k-1;                                           // j2 = k-1;
	    else                                                  // else 
	     j2 = pl-r;                                           // j2 = pl-r;
                                                                  // end

	    for (j=j1; j <= j2; j++)                              // for j = j1:j2
	      {
		a(s2,j) = (a(s1,j) - a(s1,j-1))/ndu(pk+1,rk+j);   // a(s2+1,j+1) = (a(s1+1,j+1) - a(s1+1,j))/ndu(pk+2,rk+j+1);
		d += a(s2,j)*ndu(rk+j,pk);                        // d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1);
	      }                                                   // end

	    if (r <= pk)                                          // if (r <= pk)
	      {
		a(s2,k) = -a(s1,k-1)/ndu(pk+1,r);                 // a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1);
		d += a(s2,k)*ndu(r,pk);                           // d = d + a(s2+1,k+1)*ndu(r+1,pk+1);
	      }	                                                  // end
	    
	    ders(k,r) = d;                                        // ders(k+1,r+1) = d;
	    j = s1;                                               // j = s1;
	    s1 = s2;                                              // s1 = s2;
	    s2 = j;                                               // s2 = j;
	  }                                                       // end
    }                                                             // end

  r = pl;                                                         // r = pl;
  for (k=1; k <= nders; k++)                                      // for k = 1:nders
    {
      for (j=0; j<=pl; j++)                                       // for j = 0:pl
	ders(k,j) = ders(k,j)*r;                                  // ders(k+1,j+1) = ders(k+1,j+1)*r;
                                                                  // end
      r = r*(pl-k);                                               // r = r*(pl-k);
    }                                                             // end
  
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
