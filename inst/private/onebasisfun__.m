function N = onebasisfun__ (u, p, U)

%  __ONEBASISFUN__: Undocumented internal function
%
%   Copyright (C) 2009 Carlo de Falco
%   Copyright (C) 2012 Rafael Vazquez
%   This software comes with ABSOLUTELY NO WARRANTY; see the file
%   COPYING for details.  This is free software, and you are welcome
%   to distribute it under the conditions laid out in COPYING.

  N = zeros (size (u));

  for ii = 1:numel (u)
    if (~ any (U <= u(ii))) || (~ any (U > u(ii)))
      continue;
    elseif (p == 0)
      N(ii) = 1;
      continue;
    elseif (p == 1)
      if (u(ii) < U(2))
        N(ii) = (u(ii) - U(1)) / (U(2) - U(1));
        continue;
      else
        N(ii) = (U(end) - u(ii)) / (U(3) - U(2));
        continue;
      end

    elseif (p == 2)
      if (u(ii) < U(2))
        N(ii) = (u(ii) - U(1))^2 / ((U(3) - U(1)) * (U(2) - U(1)));
        continue;
      elseif (u(ii) > U(3))
        N(ii) = (U(4) - u(ii))^2 / ((U(4) - U(2)) * (U(4) - U(3)));
        continue;
      else
        ld = U(3) - U(1); dd = U(4) - U(2);
        if (ld ~= 0)
          N(ii) = N(ii) + (u(ii) - U(1))*(U(3) - u(ii)) / ((U(3) - U(2)) * ld);
        end
        if (dd ~= 0)
          N(ii) = N(ii) + (U(4) - u(ii))*(u(ii) - U(2)) / ((U(3) - U(2)) * dd);
        end
      end

    else
      ln = u(ii) - U(1);
      ld = U(end-1) - U(1);
      if (ld ~= 0)
        N(ii) = N(ii) + ln * onebasisfun__ (u(ii), p-1, U(1:end-1))/ ld; 
      end

      dn = U(end) - u(ii);
      dd = U(end) - U(2);
      if (dd ~= 0)
        N(ii) = N(ii) + dn * onebasisfun__ (u(ii), p-1, U(2:end))/ dd;
      end
    end
  end
  
end
