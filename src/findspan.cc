/* Copyright (C) 2009 Carlo de Falco
  
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

#include <octave/oct.h>
#include "low_level_functions.h"

DEFUN_DLD(findspan, args, nargout,"\
FINDSPAN: Find the span of a B-Spline knot vector at a parametric point\n\
\n\
\n\
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
  Algorithm A2.1 from 'The NURBS BOOK' pg68\n\
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

