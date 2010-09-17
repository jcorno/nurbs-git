% KNTBRKDEGMULT: Construct an open knot vector by giving the sequence of
%                knots, the degree and the multiplicity.
%
%   knots = kntbrkdegreg (breaks, degree)
%   knots = kntbrkdegreg (breaks, degree, mult)
%
% INPUT:
%
%     breaks:  sequence of knots.
%     degree:  polynomial degree of the splines associated to the knot vector.
%     mult:    multiplicity of the knots.
%
% OUTPUT:
%
%     knots:  knot vector.
%
% If MULT has as many entries as BREAKS, or as the number of interior
%   knots, a different multiplicity will be assigned to each knot. If
%   MULT is not present, it will be taken equal to 1.
%
% Copyright (C) 2010 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function knots = kntbrkdegmult (breaks, degree, mult)

if (nargin == 2)
  mult = 1;
end

if (numel (breaks) < 2)
  error ('kntbrkdegmult: the knots sequence should contain at least two points')
end

if (numel (mult) == 1)
  mults = [degree+1, mult(ones (1, numel (breaks) - 2)), degree+1];
elseif (numel (mult) == numel (breaks))
  mults = [degree+1 mult(2:end-1) degree+1];
elseif (numel (mult) == numel (breaks) - 2)
  mults = [degree+1 mult degree+1];
else
  error('kntbrkdegmult: the length of mult should be equal to one or the number of knots')
end

if (any (mults > degree+1))
  warning ('kntbrkdegmult: some knots have higher multiplicity than the degree+1')
end

breaks = sort (breaks);

lm = numel (mults);
sm = sum (mults);

mm = zeros (1,sm);
mm (cumsum ([1 reshape(mults (1:end-1), 1, lm-1)])) = ones (1,lm);
knots = breaks (cumsum (mm));

end
