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

function sv = findspanc (uv, U)                

% FINDSPANC:  Find the span of a B-Spline knot vector at a parametric point
%
% Calling Sequence:
% 
%   s = findspanc(n,p,u,U)
% 
%  INPUT:
% 
%    u - parametric point
%    U - knot sequence
% 
%    U(1) <= u <= U(end)
%  RETURN:
% 
%    s - knot span
% 
% NOTE: the difference between findspan and findspanc is that
%        the latter works for non-open knot vectors


sv = lookup (U, uv, "lr") - 1;

%!test
%!  n = 3; 
%!  U = [0 0 0 1/2 1 1 1]; 
%!  p = 2; 
%!  u = linspace(0, 1, 10)(2:end-1);  
%!  s = findspanc (p, u, U); 
%!  assert (s, [2*ones(1, 4) 3*ones(1, 4)]);

%!test
%!  p = 2; m = 7; n = m - p - 1;
%!  U = [zeros(1,p)  linspace(0,1,m+1-2*p) ones(1,p)];
%!  u = [ 0.11880   0.55118   0.93141   0.40068   0.35492 0.44392   0.88360   0.35414   0.92186   0.83085 ];
%!  s = [2   3   4   3   3   3   4   3   4   4];
%!  assert (findspanc (p, u, U), s, 1e-10);