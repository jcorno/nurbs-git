## Copyright (C) 2003 Mark Spink, 2007 Daniel Claxton, 2009 Carlo de Falco
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
## along with this program; if not, see <http://www.gnu.org/licenses/>.

function B = basisfun (iv, uv, p, U)

% BASISFUN  Basis function for B-Spline
% -------------------------------------------------------------------------
% ADAPTATION of BASISFUN from C Routine
% -------------------------------------------------------------------------
%
% Calling Sequence:
% 
%   N = basisfun(i,u,p,U)
%   
%    INPUT:
%   
%      i - knot span  ( from FindSpan() )
%      u - parametric point
%      p - spline degree
%      U - knot sequence
%   
%    OUTPUT:
%   
%      N - Basis functions vector(numel(u)x(p+1))
%   
%    Algorithm A2.2 from 'The NURBS BOOK' pg70.

B = zeros(numel(uv), p+1);
                                               
for jj = 1:numel(uv) 

  i = iv(jj) + 1; ## findspan uses 0-based numbering
  u = uv(jj);

  left = zeros(p+1,1);
  right = zeros(p+1,1);
  
  N(1) = 1;            
  for j=1:p                    
    left(j+1) = u - U(i+1-j);
    right(j+1) = U(i+j) - u; 
    saved = 0;               
    
    for r=0:j-1              
      temp = N(r+1)/(right(r+2) + left(j-r+1));
      N(r+1) = saved + right(r+2)*temp;        
      saved = left(j-r+1)*temp;                
    end                                          
    
    N(j+1) = saved;                              
  end                                              
  
  B (jj, :) = N;

end


%!test
%!  n = 3; 
%!  U = [0 0 0 1/2 1 1 1]; 
%!  p = 2; 
%!  u = linspace (0, 1, 10);  
%!  s = findspan (n, p, u, U);  
%!  Bref = [1.00000   0.00000   0.00000
%!          0.60494   0.37037   0.02469
%!          0.30864   0.59259   0.09877
%!          0.11111   0.66667   0.22222
%!          0.01235   0.59259   0.39506
%!          0.39506   0.59259   0.01235
%!          0.22222   0.66667   0.11111
%!          0.09877   0.59259   0.30864
%!          0.02469   0.37037   0.60494
%!          0.00000   0.00000   1.00000];
%!  B = basisfun (s, u, p, U);
%!  assert (B, Bref, 1e-5);
