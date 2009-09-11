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

%% usage: skl = nrbsurfderiveval (srf, [u; v], d) 
%%   OUTPUT : skl (i, j, k, l) = i-th component derived j-1,k-1 times at the
%%   l-th point.
%%
%% Adaptation of algorithm A4.4 from the NURBS book

function skl = nrbsurfderiveval (srf, uv, d) 
  
  for iu = 1:size(uv, 2);
    wders = squeeze  (surfderiveval (srf.number(1)-1, srf.order(1)-1,  ...
				     srf.knots{1}, srf.number(2)-1,  ...
				     srf.order(2)-1, srf.knots{2},  ...
				     squeeze (srf.coefs(4, :, :)), uv(1,iu), ...
				     uv(2,iu), d));      

    for idim = 1:3
      Aders = squeeze  (surfderiveval (srf.number(1)-1, srf.order(1)-1,  ...
			srf.knots{1}, srf.number(2)-1,  ...
			srf.order(2)-1, srf.knots{2},  ...
			squeeze (srf.coefs(idim, :, :)), uv(1,iu),
			uv(2,iu), d));      
      
      for k=0:d
	for l=0:d-k
	  v = Aders(k+1, l+1);
	  for j=1:l
	    v -= bincoeff(l,j)*wders(1,j+1)*skl(idim, k+1, l-j+1,iu);
	  endfor
	  for i=1:k
	    v -= bincoeff(k+1,i+1)*wders(i+1,1)*skl(idim, k-i+1, l+1, iu);
	    v2 =0;
	    for j=1:l
	      v2 += bincoeff(l+1,j+1)*wders(i+1);
	    endfor
	    v -= bincoeff(k,i)*v2;
	  endfor
	  skl(idim, k+1, l+1, iu) = v/wders(1,1);
	endfor
      endfor
    endfor
    
  endfor

endfunction

%!test
%! k = [0 0  1 1];
%! c = [0 1];
%! [coef(2,:,:), coef(1,:,:)] = meshgrid (c, c);
%! coef(3,:,:) = coef(1,:,:);
%! srf = nrbmak (coef, {k, k});
%! uv  = linspace(0,1,11)([1 1],:);
%! skl = nrbsurfderiveval (srf, uv, 1);
%! assert (squeeze (skl (2,:,:,3)), [.2 1; 0 0])

%!shared uv, P, dPdx, d2Pdx2, c1, c2
%!test
%! uv = [sort(rand(1,10)); (rand(1,10))];
%! c1 = nrbmak([0 1/2 1; 0 1 0],[0 0 0 1 1 1]);
%! c1 = nrbtform (c1, vecrotx (pi/2));
%! c2  = nrbtform(c1, vectrans([0 1 0]));
%! srf = nrbdegelev (nrbruled (c1, c2), [0, 1]);
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
%!assert(P(3,:), 2*(P(1,:)-P(1,:).^2),100*eps)
%!assert(dPdx(3,:), 2-4*P(1,:), 100*eps)
%!assert(d2Pdx2(3,:), -4+0*P(1,:), 100*eps)

%!shared dPdu, d2Pdu2, P
%!test
%! uv  = [sort(rand(1,10)); (rand(1,10))];
%! c1 = nrbmak([0 1/2 1; 0.1 1.6 1.1; 0 0 0],[0 0 0 1 1 1]);
%! c2 = nrbmak([0 1/2 1; 0.1 1.6 1.1; 1 1 1],[0 0 0 1 1 1]);
%! srf = nrbdegelev (nrbruled (c1, c2), [0, 1]);
%! skl = nrbsurfderiveval (srf, uv, 2);
%! P = squeeze(skl(:,1,1,:));
%! dPdu = squeeze(skl(:,2,1,:));
%! dPdv = squeeze(skl(:,1,2,:));
%! d2Pdu2 = squeeze(skl(:,3,1,:));
%!assert(dPdu(2,:), 3-4*P(1,:),100*eps)
%!assert(d2Pdu2(2,:), -4+0*P(1,:),100*eps)