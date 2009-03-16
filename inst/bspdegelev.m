## Copyright (C) 2003 Mark Spink, 2007 Daniel Claxton
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

function N = basisfun(i,u,p,U)                

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
%      N - Basis functions vector[p+1]
%   
%    Algorithm A2.2 from 'The NURBS BOOK' pg70.
                                                
                                                  %   void basisfun(int i, double u, int p, double *U, double *N) {
                                                  %   int j,r;
                                                  %   double saved, temp;
i = i + 1;
                                                  %   // work space
left = zeros(p+1,1);                              %   double *left  = (double*) mxMalloc((p+1)*sizeof(double));
right = zeros(p+1,1);                             %   double *right = (double*) mxMalloc((p+1)*sizeof(double));
                                               
N(1) = 1;                                         %   N[0] = 1.0;
for j=1:p                                         %   for (j = 1; j <= p; j++) {
    left(j+1) = u - U(i+1-j);                     %   left[j]  = u - U[i+1-j];
    right(j+1) = U(i+j) - u;                      %   right[j] = U[i+j] - u;
    saved = 0;                                    %   saved = 0.0;

    for r=0:j-1                                   %   for (r = 0; r < j; r++) {
        temp = N(r+1)/(right(r+2) + left(j-r+1)); %   temp = N[r] / (right[r+1] + left[j-r]);
        N(r+1) = saved + right(r+2)*temp;         %   N[r] = saved + right[r+1] * temp;
        saved = left(j-r+1)*temp;                 %   saved = left[j-r] * temp;
    end                                           %   }

    N(j+1) = saved;                               %   N[j] = saved;
end                                               %   }
  
                                                  %   mxFree(left);
                                                  %   mxFree(right);
                                                  %   }
