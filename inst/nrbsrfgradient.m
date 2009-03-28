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
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{dzdx}, @var{dzdy}]=} nrbsrfgradient (@var{nrb}, @var{nrbder}, @var{u}, @var{v})
## Compute the gradient of a NURBS surface.
## @seealso{nrbderiv}
## @end deftypefn

## Author: Carlo de Falco <carlo@guglielmo.local>
## Created: 2009-03-17

function [dzdx, dzdy]  = nrbsrfgradient (nrb, nrbder, u, v)

  if ((nargin != 4) || (nargout>2))
    print_usage();
  endif
    
  [np, dp] = nrbdeval (nrb, nrbder, {u, v});

  dxdu = squeeze (dp{1}(1,:,:));
  dydu = squeeze (dp{1}(2,:,:));
  dzdu = squeeze (dp{1}(3,:,:));
  dxdv = squeeze (dp{2}(1,:,:));
  dydv = squeeze (dp{2}(2,:,:));
  dzdv = squeeze (dp{2}(3,:,:));

  detjac = dxdu.*dydv - dydu.*dxdv;
  
  dzdx = ( dydv .* dzdu - dydu .*dzdv)./detjac;
  dzdy = (-dxdv .* dzdu + dxdu .*dzdv)./detjac;
  
endfunction


%!shared nrb, cntl, k, rec
%!
%!test
%! p = 5; n = 5; x = y = linspace(0,3,n+1); [cntl(1,:,:), cntl(2,:,:)] = meshgrid(x,y);
%! cntl(4,:,:) = 1;
%! k{1} = k{2} = [zeros(1,p) linspace(0,1,n-p+2) ones(1,p)];
%! rec = nrbmak(cntl,k);
%! rec = nrbtform(rec, vecrotz(pi/3));
%! rec.coefs(3,:,:) = rec.coefs(1,:,:);
%! recd = nrbderiv(rec);
%! [P, w] = nrbeval(rec, {linspace(0,1,20),
%!      linspace(0,1,15)});
%! [dzdx, dzdy]  = nrbsrfgradient (rec, recd, linspace(0,1,20),
%!      linspace(0,1,15));
%! quiver(squeeze(P(1,:,:)), squeeze(P(2,:,:)), dzdx, dzdy)
%! assert(dzdx, ones(20,15),  10*eps);
%! assert(dzdy, zeros(20,15), 10*eps);
%!
%!test 
%! rec = nrbtform(rec, vecrotz(pi/2));
%! recd = nrbderiv(rec);
%! [P, w] = nrbeval(rec, {linspace(0,1,20),
%!      linspace(0,1,15)});
%! [dzdx, dzdy]  = nrbsrfgradient (rec, recd, linspace(0,1,20),
%!      linspace(0,1,15));
%! quiver(squeeze(P(1,:,:)), squeeze(P(2,:,:)), dzdx, dzdy)
%! assert(dzdy, ones(20,15),  10*eps);
%! assert(dzdx, zeros(20,15), 10*eps);
%!
%!test 
%! rec = nrbtform(rec, vecrotz(-pi/4));
%! recd = nrbderiv(rec);
%! [P, w] = nrbeval(rec, {linspace(0,1,20),
%!      linspace(0,1,15)});
%! [dzdx, dzdy]  = nrbsrfgradient (rec, recd, linspace(0,1,20),
%!      linspace(0,1,15));
%! quiver(squeeze(P(1,:,:)), squeeze(P(2,:,:)), dzdx, dzdy)
%! assert(dzdy, dzdx,  10*eps);

