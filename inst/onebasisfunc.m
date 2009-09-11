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

## usage:
## N = onebasisfunc (n, u, i, p, U)

## Author: Carlo de Falco <cdf _AT_ users.sourceforge.net>
## Created: 2009-06-25

function N = onebasisfunc (n, u, i, p, U)

  m   = length (U)-1;
  n   = mod (n, m);

  spans = mod (n:n+p, m);
  splen = U(spans+2) - U(spans+1);

  N = 0;
  if (! any (i == spans))
    return;
  elseif (p == 0)
      N = 1;
      return;
  endif

  span = find (i == spans);

  if ((span != 1)  ||  (u > U(spans(span)+1)))

    ld = sum (splen(1:p));
    if (ld != 0)
      ln = sum (splen(1:span-1)) + max (u - U(spans(span)+1), 0);    
      Nl = onebasisfun (n, u, i, p - 1, U);
      N += Nl * ln/ld;
    endif
  endif

  if ((span != p+1) || (u < U(spans(span)+2)))    
    rd = sum (splen(end-p+1:end));
    if (rd != 0)
      rn = sum (splen(span+1:end)) + max (U(spans(span)+2)-u, 0);    
      Nr = onebasisfun (n+1, u, i, p - 1, U);
      N += Nr * rn/rd;
    endif
  endif
  
endfunction

%!demo
%! p = 1;
%! U = sort([linspace(0, 1, 11) .5 .5 .5 ]);
%! uv= linspace (0, 1, 201);
%! for num = 0 %:length(U)
%!   ii = 1;
%!   for u = uv
%!     s = findspanc (u, U);
%!     N(ii++) = onebasisfunc (num, u, s, p, U);
%!   endfor
%!   plot (uv, N, [num2str(mod(num,6)) "-"]);
%!   hold on;
%! endfor
%! 
%! hold off
%!