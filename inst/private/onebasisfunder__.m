function N = onebasisfunder__ (u, p, U)

%  __ONEBASISFUNDER__: Undocumented internal function
%
%   Copyright (C) 2012 Rafael Vazquez
%   This software comes with ABSOLUTELY NO WARRANTY; see the file
%   COPYING for details.  This is free software, and you are welcome
%   to distribute it under the conditions laid out in COPYING.

  N = zeros (size (u));

  for ii = 1:numel (u)
    if (~ any (U <= u(ii))) || (~ any (U > u(ii)))
      continue;
    elseif (p == 0)
      N(ii) = 0;
      continue;
    else
      ld = U(end-1) - U(1);
      if (ld ~= 0)
        N(ii) = N(ii) + p * onebasisfun__ (u(ii), p-1, U(1:end-1))/ ld;
      end

      dd = U(end) - U(2);
      if (dd ~= 0)
        N(ii) = N(ii) - p * onebasisfun__ (u(ii), p-1, U(2:end))/ dd;
      end
    end
  end
  
end
