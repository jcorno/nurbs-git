function sv = findspan (n, p, uv, U)                
% FINDSPAN  Find the span of a B-Spline knot vector at a parametric point
% -------------------------------------------------------------------------
% ADAPTATION of FINDSPAN from C
% -------------------------------------------------------------------------
%
% Calling Sequence:
% 
%   s = findspan(n,p,u,U)
% 
%  INPUT:
% 
%    n - number of control points - 1
%    p - spline degree
%    u - parametric point
%    U - knot sequence
% 
%    U(1) <= u <= U(end)
%  RETURN:
% 
%    s - knot span
% 
%  Algorithm A2.1 from 'The NURBS BOOK' pg68
                    
if ((nargin ~= 4) || any(uv<U(1)) || any(uv>U(end))) 
  print_usage ()
end

sv = zeros(1, numel(uv));

for ii = 1:numel(uv)
  u = uv(ii);
  if (u>=U(n+2));
    s=n; 
  else
    low = p;                                        
    high = n + 1;                                   
    mid = floor((low + high) / 2);                  
    while (u < U(mid+1) || u >= U(mid+2))           
      if (u < U(mid+1))                           
	high = mid;                             
      else                                        
	low = mid;                              
      end 
      mid = floor((low + high) / 2);              
    end                                             
    s = mid;                                        
  end
  sv(ii) = s;
end

%!test
%!  n = 3; 
%!  U = [0 0 0 1/2 1 1 1]; 
%!  p = 2; 
%!  u = linspace(0, 1, 10);  
%!  s = findspan (n, p, u, U); 
%!  assert (s, [2*ones(1, 5) 3*ones(1, 5)]);

