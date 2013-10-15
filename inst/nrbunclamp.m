function ucrv = nrbunclamp (crv, k)

% NRBUNCLAMP: Compute the knot vector and control points of the unclamped curve
%
% Calling Sequence:
% 
%   ucrv = nrbrunclamp (crv, k)
% 
% INPUT:
% 
%   crv	: NURBS curve, see nrbmak.
% 
%   k   : continuity for the unclamping (from 0 up to p-1)
%
% OUTPUT:
% 
%   ucrv: NURBS curve with unclamped knot vector, see nrbmak
% 
% Description:
% 
%   Unclamps a curve, removing the open knot vector. Computes the new
%     knot vector and control points of the unclamped curve.
%
%    Adapted from Algorithm A12.1 from 'The NURBS BOOK' pg577.
% 
%    Copyright (C) 2013 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

  p  = crv.order - 1;
  U  = crv.knots;
  n  = crv.number;
  Pw = crv.coefs;
  m = n + p + 1;

  if (k >= p)
    warning ('Taking the maximum k allowed, degree - 1')
    k = p - 1;
  end

  % Unclamp at left end
  for ii=0:k
    U(k-ii+1) = U(k-ii+2) - (U(n+1-ii) - U(n-ii));
  end

  for ii = p-k-1:p-2
    for jj = ii:-1:0
      alpha = (U(p+1) - U(p+jj-ii)) / (U(p+jj+2) - U(p+jj-ii));
      Pw(:,jj+1) = (Pw(:,jj+1) - alpha*Pw(:,jj+2))/(1-alpha);
    end
  end


  % Unclamp at right end
  for ii=0:k
    U(m-k+ii) = U(m-k+ii-1) + U(p+ii+1+1) - U(p+ii+1);
  end

  for ii = p-k-1:p-2
    for jj = ii:-1:0
      alpha = (U(n+1)-U(n-jj))/(U(n+2-jj+ii)-U(n-jj));
      Pw(:,n-jj) = (Pw(:,n-jj) - (1-alpha)*Pw(:,n-jj-1))/alpha;
    end
  end

ucrv = nrbmak (Pw, U);
  
%!demo
%! crv = nrbcirc (1,[],0,2*pi/3);
%! crv = nrbdegelev (crv, 2);
%! figure
%! nrbctrlplot (crv); hold on
%! nrbctrlplot (nrbtform (nrbunclamp (crv, 1), vectrans([-0.4, -0.4])));
%! nrbctrlplot (nrbtform (nrbunclamp (crv, 2), vectrans([-0.8, -0.8])));
%! nrbctrlplot (nrbtform (nrbunclamp (crv, 3), vectrans([-1.6, -1.6])));
%! title ('Original curve and unclamped versions')

%!test
%! crv = nrbdegelev (nrbtestcrv,2);
%! x = linspace (0, 1, 100);
%! F = nrbeval (crv, x);
%! ucrv = nrbunclamp (crv, 0);
%! assert (F, nrbeval(ucrv, x));
%! ucrv = nrbunclamp (crv, 1);
%! assert (F, nrbeval(ucrv, x), 1e-14);
%! ucrv = nrbunclamp (crv, 2);
%! assert (F, nrbeval(ucrv, x), 1e-14);
%! ucrv = nrbunclamp (crv, 3);
%! assert (F, nrbeval(ucrv, x), 1e-14);
