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

function varargout = nrbbasisfunder (points, nrb)

% NRBBASISFUNDER:  NURBS basis functions derivatives
%
% Calling Sequence:
% 
%   Bu          = nrbbasisfunder (u, crv)
%   [Bu, N]     = nrbbasisfunder (u, crv)
%   [Bu, Bv]    = nrbbasisfunder ({u, v}, srf)
%   [Bu, Bv, N] = nrbbasisfunder ({u, v}, srf)
%
%    INPUT:
%   
%      u   - parametric points along u direction
%      v   - parametric points along v direction
%      crv - NURBS curve
%      srf - NURBS surface
%   
%    OUTPUT:
%   
%      Bu - Basis functions derivatives WRT direction u
%           size(Bu)=[numel(u),(p+1)] for curves
%           or [numel(u)*numel(v), (p+1)*(q+1)] for surfaces
%
%      Bv - Basis functions derivatives WRT direction v
%           size(Bv)=[numel(v),(p+1)] for curves
%           or [numel(u)*numel(v), (p+1)*(q+1)] for surfaces
%
%      N - Indices of the basis functions that are nonvanishing at each
%          point. size(N) == size(B)
%   

  
  if (   (nargin<2) 
      || (nargout>3)
      || (~isstruct(nrb))
      || (iscell(points) && ~iscell(nrb.knots))
      || (~iscell(points) && iscell(nrb.knots) && (size(points,1)~=2))
      || (~iscell(nrb.knots) && (nargout>2))
      )
    print_usage();
  end
                            
  if (~iscell(nrb.knots))    %% NURBS curve
    
    [varargout{1}, varargout{2}] = __nrb_crv_basisfun_der__ (points, nrb);

  else                       %% NURBS surface

    [varargout{1}, varargout{2}, varargout{3}] = __nrb_srf_basisfun_der__ (points, nrb);

  end

  function [Bu, nbfu] = __nrb_crv_basisfun_der__ (points, nrb);
    n    = size (nrb.coefs, 2) -1;
    p    = nrb.order -1;
    u    = points;
    U    = nrb.knots;
    w    = nrb.coefs(4,:);
    
    spu  =  findspan (n, p, u, U);
    nbfu =  numbasisfun (spu, u, p, U);
    
    N     = basisfun (spu, u, p, U);

    Nprime = basisfunder (spu, p, u, U, 1);
    Nprime = squeeze(Nprime(:,2,:));

    
    [Dpc, Dpk]  = bspderiv (p, w, U);
    D           = bspeval  (p, w, U, u);
    Dprime      = bspeval  (p-1, Dpc, Dpk, u);
    

    Bu1   = bsxfun (@(np, d) np/d , Nprime.', D);
    Bu2   = bsxfun (@(n, dp)  n*dp, N.', Dprime./D.^2);
    Bu    = w(nbfu+1) .* (Bu1 - Bu2).';

  end

  function [Bu, Bv, N] = __nrb_srf_basisfun_der__ (points, nrb);

    m    = size (nrb.coefs, 2) -1;
    n    = size (nrb.coefs, 3) -1;
    
    p    = nrb.order(1) -1;
    q    = nrb.order(2) -1;

    if (iscell(points))
      [v, u] = meshgrid(points{2}, points{1});
      u = reshape(u, 1, []); v = reshape(v, 1, []);
    else
      u = points(1,:,:)(:).';
      v = points(2,:,:)(:).';
    end

    npt = length(u);

    U    = nrb.knots{1};
    V    = nrb.knots{2};
    
    w    = squeeze(nrb.coefs(4,:,:));

    spu  =  findspan (m, p, u, U); 
    Ik   =  numbasisfun (spu, u, p, U);

    spv  =  findspan (n, q, v, V);
    Jk   =  numbasisfun (spv, v, q, V);
    
    NuIkuk = basisfun (spu, u, p, U);
    NvJkvk = basisfun (spv, v, q, V);

    NuIkukprime = basisfunder (spu, p, u, U, 1);
    NuIkukprime = squeeze(NuIkukprime(:,2,:));

    NvJkvkprime = basisfunder (spv, q, v, V, 1);
    NvJkvkprime = squeeze(NvJkvkprime(:,2,:));

    
    for k=1:npt
      [Ika, Jkb] = meshgrid(Ik(k, :), Jk(k, :)); 
      
      N(k, :)    = sub2ind([m+1, n+1], Ika(:)+1, Jkb(:)+1);
      wIkaJkb(1:p+1, 1:q+1) = reshape (w(N(k, :)), p+1, q+1); 

      Num    = (NuIkuk(k, :).' * NvJkvk(k, :)) .* wIkaJkb;
      Num_du = (NuIkukprime(k, :).' * NvJkvk(k, :)) .* wIkaJkb;
      Num_dv = (NuIkuk(k, :).' * NvJkvkprime(k, :)) .* wIkaJkb;
      Denom  = sum(sum(Num));
      Denom_du = sum(sum(Num_du));
      Denom_dv = sum(sum(Num_dv));

      Bu(k, :) = (Num_du/Denom - Denom_du.*Num/Denom.^2)(:).';
      Bv(k, :) = (Num_dv/Denom - Denom_dv.*Num/Denom.^2)(:).';
    end
    
  end
  
  

%!test
%! U = [0 0 0 0 1 1 1 1];
%! x = [0 1/3 2/3 1] ;
%! y = [0 0 0 0];
%! w = rand(1,4);
%! nrb = nrbmak ([x;y;y;w], U);
%! u = linspace(0, 1, 30);
%! [Bu, id] = nrbbasisfunder (u, nrb);
%! #plot(u, Bu)
%! assert (sum(Bu, 2), zeros(numel(u), 1), 1e-10), 

%!test
%! U = [0 0 0 0 1/2 1 1 1 1];
%! x = [0 1/4 1/2 3/4 1] ;
%! y = [0 0 0 0 0];
%! w = rand(1,5);
%! nrb = nrbmak ([x;y;y;w], U);
%! u = linspace(0, 1, 300); 
%! [Bu, id] = nrbbasisfunder (u, nrb); 
%! assert (sum(Bu, 2), zeros(numel(u), 1), 1e-10)

%!test
%! p = 2;
%! q = 3;
%! mcp = 2; ncp = 3;
%! knots = {[zeros(1,p), linspace(0,1,mcp-p+2), ones(1,p)], [zeros(1,q), linspace(0,1,ncp-q+2), ones(1,q)]};
%! Lx  = 1; Ly  = 1;
%! [cntl(2,:,:), cntl(1,:,:)] = meshgrid(linspace(0, Ly, ncp+1), linspace(0, Lx, mcp+1) );
%! cntl(4,:,:) = 1:numel(cntl(1,:,:));
%! nrb = nrbmak(cntl, knots);
%! [u(1,:,:), u(2,:,:)] = meshgrid(rand (1, 20), rand (1, 20));
%! tic(); [Bu, Bv, N] = nrbbasisfunder (u, nrb); toc()
%! #plot3(squeeze(u(1,:,:)), squeeze(u(2,:,:)), reshape(Bu(:,10), 20, 20),'o')
%! assert (sum (Bu, 2), zeros(20^2, 1), 1e-10)

