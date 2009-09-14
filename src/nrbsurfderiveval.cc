/*
%% Copyright (C) 2009 Carlo de Falco
%% 
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%% 
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%% 
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <octave/oct.h>
#include <octave/oct-map.h>
#include "low_level_functions.h"

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
   

DEFUN_DLD(nrbsurfderiveval, args, nargout,"\
\n usage: skl = nrbsurfderiveval (srf, [u; v], d) \
\n   OUTPUT : skl (i, j, k, l) = i-th component derived j-1,k-1 times at the\
\n   l-th point.\
\n\
\n Adaptation of algorithm A4.4 from the NURBS book\n")
{
  //function skl = nrbsurfderiveval (srf, uv, d) 
  
  octave_value_list retval;

  Octave_map srf = args(0).map_value();
  Matrix uv = args(1).matrix_value ();
  octave_idx_type d = args(2).idx_type_value ();

  if (! error_state)
    {
      Array<octave_idx_type> idxta (4, 0);
      dim_vector idxa; idxa.resize (4);
      idxa(0) = 3; idxa(1) = d+1; 
      idxa(2) = d+1; idxa(3) = uv.columns (); 
      NDArray skl (idxa, 0.0);

      for (octave_idx_type iu(0); iu<uv.cols (); iu++)
	{
	  octave_idx_type n = octave_idx_type ((srf.contents("number")(0).row_vector_value())(0) - 1);
	  octave_idx_type m = octave_idx_type ((srf.contents("number")(0).row_vector_value())(1) - 1);
	  octave_idx_type p = octave_idx_type ((srf.contents("order")(0).row_vector_value())(0) - 1);
	  octave_idx_type q = octave_idx_type ((srf.contents("order")(0).row_vector_value())(1) - 1);
	  Cell knots = srf.contents("knots")(0).cell_value();
	  RowVector knotsu = knots.elem (0).row_vector_value ();
	  RowVector knotsv = knots.elem (1).row_vector_value ();
	  NDArray coefs  = srf.contents("coefs")(0).array_value();
	  Array<idx_vector> idx(3, idx_vector(':'));
	  idx (0) = idx_vector (3);
	  Matrix weights (NDArray (coefs.index (idx).squeeze ()).matrix_value ());
	  
	  Matrix wders;
	  surfderiveval (n, p, knotsu, m, q, knotsv, weights, uv(0,iu), uv(1,iu), d, wders);      
	  
	  for (octave_idx_type idim (0); idim<=2; idim++)
	    {

	      Matrix Aders;
	      idx (0) = idx_vector (idim);
	      Matrix P (NDArray (coefs.index (idx).squeeze ()).matrix_value ());
	      //std::cout << P << "\n\n\n";
	      surfderiveval (n, p, knotsu, m, q, knotsv, P, uv(0,iu), uv(1,iu), d, Aders);;      
	      
	      for (octave_idx_type k(0); k<=d; k++)
		{
		  for (octave_idx_type l(0); l<=d-k; l++)
		    {
		      double v = Aders(k, l);
		      for (octave_idx_type j(1); j<=l; j++)
			{
			  idxta(0) = idim; idxta(1) = k; idxta(2) = l; idxta(3) = iu;
			  v -= bincoeff(l,j) * wders(0,j) * skl(idxta);
			}
		      for (octave_idx_type i(1); i<=k; i++)
			{
			  idxta(0) = idim; idxta(1) = k-i; idxta(2) = l; idxta(3) = iu;
			  v -= bincoeff(k,i) * wders(i,0) * skl(idxta);
			  double v2 = 0.0;
			  for (octave_idx_type j(1);j<=l;j++)
			    {
			      v2 += bincoeff(l,j) * wders(i);
			    }
			  v -= bincoeff(k,i) * v2;
			}
		      idxta(0) = idim; idxta(1) = k; idxta(2) = l; idxta(3) = iu;
		      skl(idxta) = v/wders(0,0);
		    }
		}
	    }
	  
	} 
      retval(0) = octave_value (skl);
    }
  return retval;
}
/*
%!test
%! k = [0 0  1 1];
%! c = [0 1];
%! [coef(2,:,:), coef(1,:,:)] = meshgrid (c, c);
%! coef(3,:,:) = coef(1,:,:);
%! srf = nrbmak (coef, {k, k});
%! [u, v] = meshgrid (linspace(0,1,11));
%! uv = [u(:)';v(:)'];
%! skl = nrbsurfderiveval (srf, uv, 0);
%! assert (squeeze (skl (1:2,1,1,:)), nrbeval (srf, uv)(1:2,:), 1e3*eps)


%!test
%! k = [0 0  1 1];
%! c = [0 1];
%! [coef(2,:,:), coef(1,:,:)] = meshgrid (c, c);
%! coef(3,:,:) = coef(1,:,:);
%! srf = nrbmak (coef, {k, k});
%! srf = nrbkntins (srf, {[], rand(2,1)});
%! [u, v] = meshgrid (linspace(0,1,11));
%! uv = [u(:)';v(:)'];
%! skl = nrbsurfderiveval (srf, uv, 0);
%! assert (squeeze (skl (1:2,1,1,:)), nrbeval (srf, uv)(1:2,:), 1e3*eps)

%!shared srf, uv
%!test 
%! k = [0 0 0 1 1 1];
%! c = [0 1/2 1];
%! [coef(1,:,:), coef(2,:,:)] = meshgrid (c, c);
%! coef(3,:,:) = coef(1,:,:);
%! srf = nrbmak (coef, {k, k});
%! ders= nrbderiv (srf);
%! [u, v] = meshgrid (linspace(0,1,11));
%! uv = [u(:)';v(:)'];
%! skl = nrbsurfderiveval (srf, uv, 1);
%! [fun, der] = nrbdeval (srf, ders, uv);
%! assert (squeeze (skl (1:2,1,1,:)), fun(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,2,1,:)), der{1}(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,1,2,:)), der{2}(1:2,:), 1e3*eps)
%!
%!test 
%! srf = nrbdegelev (srf, [3, 1]);
%! ders= nrbderiv (srf);
%! [fun, der] = nrbdeval (srf, ders, uv);
%! skl = nrbsurfderiveval (srf, uv, 1);
%! assert (squeeze (skl (1:2,1,1,:)), fun(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,2,1,:)), der{1}(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,1,2,:)), der{2}(1:2,:), 1e3*eps)

%!shared uv
%!test 
%! k = [0 0 0 1 1 1];
%! c = [0 1/2 1];
%! [coef(2,:,:), coef(1,:,:)] = meshgrid (c, c);
%! coef(3,:,:) = coef(1,:,:);
%! srf = nrbmak (coef, {k, k});
%! ders= nrbderiv (srf);
%! [u, v] = meshgrid (linspace(0,1,11));
%! uv = [u(:)';v(:)'];
%! skl = nrbsurfderiveval (srf, uv, 1);
%! [fun, der] = nrbdeval (srf, ders, uv);
%! assert (squeeze (skl (1:2,1,1,:)), fun(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,2,1,:)), der{1}(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,1,2,:)), der{2}(1:2,:), 1e3*eps)
%!
%!test 
%! p = q = 3;
%! mcp = 5; ncp = 5;
%! Lx  = Ly  = 10*rand(1);
%! srf = nrbdegelev (nrb4surf ([0 0], [Lx, 0], [0 Ly], [Lx Ly]), [p-1, q-1]);
%! %%srf = nrbkntins (srf, {linspace(0,1,mcp-p+2)(2:end-1), linspace(0,1,ncp-q+2)(2:end-1)});
%! %%srf.coefs = permute (srf.coefs, [1 3 2]);
%! ders= nrbderiv (srf);
%! [fun, der] = nrbdeval (srf, ders, uv);
%! skl = nrbsurfderiveval (srf, uv, 1);
%! assert (squeeze (skl (1:2,1,1,:)), fun(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,2,1,:)), der{1}(1:2,:), 1e3*eps)
%! assert (squeeze (skl (1:2,1,2,:)), der{2}(1:2,:), 1e3*eps)

%!shared srf, uv, P, dPdx, d2Pdx2, c1, c2
%!test
%! [u, v] = meshgrid (linspace(0,1,10));
%! uv = [u(:)';v(:)'];
%! c1 = nrbmak([0 1/2 1; 0 1 0],[0 0 0 1 1 1]);
%! c1 = nrbtform (c1, vecrotx (pi/2));
%! c2  = nrbtform(c1, vectrans([0 1 0]));
%! srf = nrbdegelev (nrbruled (c1, c2), [3, 1]);
%! skl = nrbsurfderiveval (srf, uv, 2);
%! P = squeeze(skl(:,1,1,:));
%! dPdx = squeeze(skl(:,2,1,:));
%! d2Pdx2 = squeeze(skl(:,3,1,:));
%!assert(P(3,:), 2*(P(1,:)-P(1,:).^2),100*eps)
%!assert(dPdx(3,:), 2-4*P(1,:), 100*eps)
%!assert(d2Pdx2(3,:), -4+0*P(1,:), 100*eps)
%! srf = nrbdegelev (nrbruled (c1, c2), [5, 6]);
%! skl = nrbsurfderiveval (srf, uv, 2);
%! P = squeeze(skl(:,1,1,:));
%! dPdx = squeeze(skl(:,2,1,:));
%! d2Pdx2 = squeeze(skl(:,3,1,:));
%! assert (squeeze (skl (1:2,1,1,:)), nrbeval (srf, uv)(1:2,:), 1e3*eps)
%!assert(P(3,:), 2*(P(1,:)-P(1,:).^2),100*eps)
%!assert(dPdx(3,:), 2-4*P(1,:), 100*eps)
%!assert(d2Pdx2(3,:), -4+0*P(1,:), 100*eps)
%!
%!test
%! skl = nrbsurfderiveval (srf, uv, 0);
%! assert (squeeze (skl (1:2,1,1,:)), nrbeval (srf, uv)(1:2,:), 1e3*eps)

%!shared dPdu, d2Pdu2, P, srf, uv
%!test
%! [u, v] = meshgrid (linspace(0,1,10));
%! uv = [u(:)';v(:)'];
%! c1 = nrbmak([0 1/2 1; 0.1 1.6 1.1; 0 0 0],[0 0 0 1 1 1]);
%! c2 = nrbmak([0 1/2 1; 0.1 1.6 1.1; 1 1 1],[0 0 0 1 1 1]);
%! srf = nrbdegelev (nrbruled (c1, c2), [0, 1]);
%! skl = nrbsurfderiveval (srf, uv, 2);
%! P = squeeze(skl(:,1,1,:));
%! dPdu = squeeze(skl(:,2,1,:));
%! dPdv = squeeze(skl(:,1,2,:));
%! d2Pdu2 = squeeze(skl(:,3,1,:));
%! assert (squeeze (skl (1:2,1,1,:)), nrbeval (srf, uv)(1:2,:), 1e3*eps)
%!assert(dPdu(2,:), 3-4*P(1,:),100*eps)
%!assert(d2Pdu2(2,:), -4+0*P(1,:),100*eps)
%!
%!test
%! skl = nrbsurfderiveval (srf, uv, 0);
%! assert (squeeze (skl (1:2,1,1,:)), nrbeval (srf, uv)(1:2,:), 1e3*eps)

*/
