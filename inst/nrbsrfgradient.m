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
## @deftypefn {Function File} {[@var{dzdx}, @var{dzdy}, @var{z}, @var{x}, @var{y}]=} nrbsrfgradient (@var{nrb}, @var{nrbder}, @var{u}, @var{v})
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

  dzdx = dzdy = zeros(length(u), length(v));

  dzdx(dxdu~=0) += (dzdu./dxdu)(dxdu~=0);
  dzdx(dxdv~=0) += (dzdv./dxdv)(dxdv~=0);

  dzdy(dydu~=0) += (dzdu./dydu)(dydu~=0);
  dzdy(dydv~=0) += (dzdv./dydv)(dydv~=0);

endfunction
