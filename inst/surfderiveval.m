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

%% usage: pkl = surfderivcpts (n, p, U, m, q, V, P, d) 
%%
%%  INPUT: 
%%        n+1, m+1 = number of control points
%%        p, q     = spline order
%%        U, V     = knots
%%        P        = control points
%%        u,v      = evaluation points
%%        d        = derivative order
%%  OUTPUT:
%%        skl (k+1, l+1) =  surface differentiated k
%%                          times in the u direction and l
%%                          times in the v direction
%%
%% Adaptation of algorithm A3.8 from the NURBS book

function skl = surfderiveval (n, p, U, m, q, V, P, u, v, d) 
  
  du = min (d, p);   
  for k=p+1:d
    for l=0:d-k
      skl(k+1,l+1)=0;
    endfor
  endfor

  dv = min (d, q);   
  for l=q+1:d
    for k=0:d-l
      skl(k+1,l+1)=0;
    endfor
  endfor

  uspan = findspan (n, p, u, U);
  for ip=0:p
      Nu(1:ip+1,ip+1) = basisfun (uspan, u, ip, U)';
  endfor
  
  vspan = findspan (m, q, v, V);
  for ip=0:q
      Nv(1:ip+1,ip+1) = basisfun (vspan, v, ip, V)';
  endfor

  pkl = surfderivcpts (n, p, U, m, q, V, P, d, uspan-p, uspan,  ...
		       vspan-q, vspan);

  for k = 0:du
    dd = min (d-k, dv);
    for l = 0:dd
      skl(k+1,l+1) =0;
      for i=0:q-l
	tmp = 0;
	for j = 0:p-k
	  tmp += Nu(j+1,p-k+1) * pkl(k+1,l+1,j+1,i+1);
	endfor
	skl(k+1,l+1) += Nv(i+1,q-l+1)*tmp;
      endfor
    endfor
  endfor
  
endfunction

%!shared srf
%!test
%! k = [0 0 0 1 1 1];
%! c = [0 1/2 1];
%! [coef(2,:,:), coef(1,:,:)] = meshgrid (c, c);
%! srf = nrbmak (coef, {k, k});
%! skl = surfderiveval (srf.number(1)-1, 
%!                      srf.order(1)-1, 
%!                      srf.knots{1}, 
%!                      srf.number(2)-1, 
%!                      srf.order(2)-1, 
%!                      srf.knots{2},
%!                      squeeze(srf.coefs(1,:,:)), .5, .5, 1) ;
%! assert (skl, [.5 0; 1 0])
%!test
%! srf = nrbkntins (srf, {[], rand(1,2)});
%! skl = surfderiveval (srf.number(1)-1, 
%!                      srf.order(1)-1, 
%!                      srf.knots{1},
%!                      srf.number(2)-1, 
%!                      srf.order(2)-1, 
%!                      srf.knots{2},
%!                      squeeze(srf.coefs(1,:,:)), .5, .5, 1) ;
%! assert (skl, [.5 0; 1 0], 100*eps)