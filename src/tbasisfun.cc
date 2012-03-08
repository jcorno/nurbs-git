/* Copyright (C) 2009 Carlo de Falco
   Copyright (C) 2012 Rafael Vazquez

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <octave/oct.h>
#include <iostream>

double onebasisfun__ (double u, octave_idx_type p, RowVector U)
{

  //std::cout << "u=" << u << " " << "p=" << p << " \n" << "U=" << U;
  
  double N = 0.0;
  if ((u <= U.min ()) || ( u > U.max ()))
    return (N);
  else if (p == 0)
    return (1.0);
  else if (p == 1) {
    if (u < U(1)) {
      N = (u - U(0)) / (U(1) - U(0));
      return (N);
    }
    else {
      N = (U(2) - u) / (U(2) - U(1));
      return (N);
    }
  }
  else if (p == 2) {
    double ln = u - U(0);
    double dn = U(3) - u;
    double ld = U(2) - U(0); 
    double dd = U(3) - U(1);
    if (u < U(1)) {
      N = ln*ln / (ld * (U(1) - U(0)));
      return (N);
    }
    else if (u > U(2)) {
      N = dn*dn / (dd * (U(3) - U(2)));
      return (N);
    }
    else {
      if (ld != 0)
        N = N + ln * (U(2) - u) / ((U(2) - U(1)) * ld);
      if (dd != 0)
        N = N + dn * (u - U(1)) / ((U(2) - U(1)) * dd);
      return (N);
    }
  }
 
  double ln = u - U(0);
  double ld = U(U.length () - 2) - U(0);
  if (ld != 0)
    N += ln * onebasisfun__ (u, p-1, U.extract (0, U.length () - 2))/ ld; 
    
  double dn = U(U.length () - 1) - u;
  double dd = U(U.length () - 1) - U(1);
  if (dd != 0)
    N += dn * onebasisfun__ (u, p-1, U.extract (1, U.length () - 1))/ dd;
    
  return (N);
}


double onebasisfunder__ (double u, octave_idx_type p, RowVector U)
{

  //std::cout << "u=" << u << " " << "p=" << p << " \n" << "U=" << U;
  
  double N = 0.0;
  if ((u <= U.min ()) || ( u > U.max ()))
    return (N);
  else if (p == 0)
    return (0.0);
  else {
 
  double ld = U(U.length () - 2) - U(0);
  if (ld != 0)
    N += p * onebasisfun__ (u, p-1, U.extract (0, U.length () - 2))/ ld; 
    
  double dd = U(U.length () - 1) - U(1);
  if (dd != 0)
    N -= p * onebasisfun__ (u, p-1, U.extract (1, U.length () - 1))/ dd;
    
  return (N);
  }
}

   
DEFUN_DLD(tbasisfun, args, nargout,"\
TBASISFUN: Compute a B- or T-Spline basis function, and its derivatives, from its local knot vector.\n\
\n\
 usage:\n\
\n\
 [N, Nder] = tbasisfun (u, p, U)\n\
 [N, Nder] = tbasisfun ([u; v], [p q], {U, V})\n\
 [N, Nder] = tbasisfun ([u; v; w], [p q r], {U, V, W})\n\
 \n\
 INPUT:\n\
  u or [u; v] : points in parameter space where the basis function is to be\n\
  evaluated \n\
  \n\
  U or {U, V} : local knot vector\n\
\n\
  p or [p q] : polynomial order of the basis function\n\
\n\
 OUTPUT:\n\
  N : basis function evaluated at the given parametric points\n\
  Nder : gradient of the basis function evaluated at the given points\n")

{
  
  octave_value_list retval;
  Matrix u = args(0).matrix_value ();

  RowVector N(u.cols ());
  if (! args(2).is_cell ())
    {

      double p = args(1).idx_type_value ();
      RowVector U = args(2).row_vector_value (true, true);
      assert (U.numel () == p+2);
      
      for (octave_idx_type ii=0; ii<u.numel (); ii++)
	N(ii) = onebasisfun__ (u(ii), p, U);

      if (nargout == 2) {
        RowVector Nder(u.cols ());
        for (octave_idx_type ii=0; ii<u.numel (); ii++)
          Nder(ii) = onebasisfunder__ (u(ii), p, U);
        retval(1) = Nder;
      }

    }  else {
    RowVector p = args(1).row_vector_value ();

    if (p.length() == 2) {
      Cell C = args(2).cell_value ();
      RowVector U = C(0).row_vector_value (true, true);
      RowVector V = C(1).row_vector_value (true, true);
      double Nu, Nv;
    
      for (octave_idx_type ii=0; ii<u.cols (); ii++)
        {
          Nu = onebasisfun__ (u(0, ii), octave_idx_type(p(0)), U);
          Nv = onebasisfun__ (u(1, ii), octave_idx_type(p(1)), V);
          N(ii) = Nu * Nv;
	//std::cout << "N=" << N(ii) << "\n\n\n";
        }

      if (nargout == 2) {
        Matrix Nder (2, u.cols());
        for (octave_idx_type ii=0; ii<u.cols (); ii++)
          {
            Nder(0,ii) = onebasisfunder__ (u(0, ii), octave_idx_type(p(0)), U) *
              Nv;
            Nder(1,ii) = onebasisfunder__ (u(1, ii), octave_idx_type(p(1)), V) *
	      Nu;
	//std::cout << "N=" << N(ii) << "\n\n\n";
          }
        retval(1) = Nder;
        }

    } else if (p.length() == 3) {
      Cell C = args(2).cell_value ();
      RowVector U = C(0).row_vector_value (true, true);
      RowVector V = C(1).row_vector_value (true, true);
      RowVector W = C(2).row_vector_value (true, true);
      double Nu, Nv, Nw;
    
      for (octave_idx_type ii=0; ii<u.cols (); ii++)
        {
          Nu = onebasisfun__ (u(0, ii), octave_idx_type(p(0)), U);
          Nv = onebasisfun__ (u(1, ii), octave_idx_type(p(1)), V);
          Nw = onebasisfun__ (u(2, ii), octave_idx_type(p(2)), W);
          N(ii) = Nu * Nv * Nw;
	//std::cout << "N=" << N(ii) << "\n\n\n";
        }

      if (nargout == 2) {
        Matrix Nder (3, u.cols());
        for (octave_idx_type ii=0; ii<u.cols (); ii++)
          {
            Nder(0,ii) = onebasisfunder__ (u(0, ii), octave_idx_type(p(0)), U) *
              Nv * Nw;
            Nder(1,ii) = onebasisfunder__ (u(1, ii), octave_idx_type(p(1)), V) *
              Nu * Nw;
            Nder(2,ii) = onebasisfunder__ (u(2, ii), octave_idx_type(p(2)), W) *
              Nu * Nv;
	//std::cout << "N=" << N(ii) << "\n\n\n";
          }
        retval(1) = Nder;
        }

    }
  }
  retval(0) = octave_value (N);
  return retval;
}

/*

%!demo
%! U = {[0 0 1/2 1 1], [0 0 0 1 1]};
%! p = [3, 3];
%! [X, Y] = meshgrid (linspace(0, 1, 30));
%! u = [X(:), Y(:)]';
%! N = tbasisfun (u, p, U);
%! surf (X, Y, reshape (N, size(X)))
%! title('Basis function associated to a local knot vector')
%! hold off

%!test
%! U = [0 1/2 1];
%! p = 1;
%! u = [0.3 0.4 0.6 0.7];
%! [N, Nder] = tbasisfun (u, p, U);
%! assert (N, [0.6 0.8 0.8 0.6], 1e-12);
%! assert (Nder, [2 2 -2 -2], 1e-12);

%!test
%! U = {[0 1/2 1] [0 1/2 1]};
%! p = [1 1];
%! u = [0.3 0.4 0.6 0.7; 0.3 0.4 0.6 0.7];
%! [N, Nder] = tbasisfun (u, p, U);
%! assert (N, [0.36 0.64 0.64 0.36], 1e-12);
%! assert (Nder, [1.2 1.6 -1.6 -1.2; 1.2 1.6 -1.6 -1.2], 1e-12);

%!test
%! U = {[0 1/2 1] [0 1/2 1] [0 1/2 1]};
%! p = [1 1 1];
%! u = [0.4 0.4 0.6 0.6; 0.4 0.4 0.6 0.6; 0.4 0.6 0.4 0.6];
%! [N, Nder] = tbasisfun (u, p, U);
%! assert (N, [0.512 0.512 0.512 0.512], 1e-12);
%! assert (Nder, [1.28 1.28 -1.28 -1.28; 1.28 1.28 -1.28 -1.28; 1.28 -1.28 1.28 -1.28], 1e-12);

*/
