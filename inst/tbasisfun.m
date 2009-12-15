## Copyright (C) 2009 Carlo de Falco
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## TBASISFUN: Compute a B- or T-Spline basis function from its local knot vector.
##
## usage:
##
## N = tbasisfun (u, p, U)
## N = tbasisfun ([u; v], [p q], {U, V})
## 
## INPUT:
##  u or [u; v] : points in parameter space where the basis function is to be
##  evaluated 
##  
##  U or {U, V} : local knot vector
##
## p or [p q] : polynomial order of the bais function
##
## OUTPUT:
##  N : basis function evaluated at the given parametric points 
 
## Author: Carlo de Falco <cdf _AT_ users.sourceforge.net>
## Created: 2009-06-25

function N = tbasisfun (u, p, U)
  
  if (! iscell (U))
    U = sort (U);
    assert (numel (U) == p+2)
    
    for ii=1:numel(u)
      N(ii) = onebasisfun__ (u(ii), p, U);
    endfor

  else

    for ii=1:numel(u(1,:))
      Nu(ii) = onebasisfun__ (u(1,ii), p(1), U{1});
    endfor

    for ii=1:numel(u(1,:))
      Nv(ii) = onebasisfun__ (u(2,ii), p(2), U{2});
    endfor

    N = Nu.*Nv;
  endif

endfunction

function N = onebasisfun__ (u, p, U)

  N = 0;
  if (! any (U <= u)) | (! any (U >= u))
    return;
  elseif (p == 0)
      N = 1;
      return;
  endif

 ln = u - U(1);
 ld = U(end-1) - U(1);
 if (ld != 0)
   N += ln * onebasisfun__ (u, p-1, U(1:end-1))/ ld; 
 endif

 dn = U(end) - u;
 dd = U(end) - U(2);
 if (dd != 0)
   N += dn * onebasisfun__ (u, p-1, U(2:end))/ dd;
 endif
  
endfunction

%!demo
%! U = {[0 0 1/2 1 1], [0 0 0 1 1]};
%! p = [3, 3];
%! [X, Y] = meshgrid (linspace(0, 1, 30));
%! u = [X(:), Y(:)]';
%! N = tbasisfun (u, p, U);
%! surf (X, Y, reshape (N, size(X)))
