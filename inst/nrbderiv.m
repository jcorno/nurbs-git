function varargout = nrbderiv (nurbs)
% 
% NRBDERIV2: Construct the first and second derivative representation of a
%           NURBS curve, surface or volume.
% 
% Calling Sequence:
% 
%   ders = nrbderiv (nrb);
%   [ders, ders2] = nrbderiv (nrb);
% 
% INPUT:
% 
%   nrb		: NURBS data structure, see nrbmak.
%
% OUTPUT:
% 
%   ders:  A data structure that represents the first
% 		    derivatives of a NURBS curve, surface or volume.
%   ders2: A data structure that represents the second
% 		    derivatives of a NURBS curve, surface or volume.
% 
% Description:
% 
%   The derivatives of a B-Spline are themselves a B-Spline of lower degree,
%   giving an efficient means of evaluating multiple derivatives. However,
%   although the same approach can be applied to NURBS, the situation for
%   NURBS is more complex. We have followed in this function the same idea
%   that was already used for the first derivative in the function nrbderiv.
%   The second derivative data structure can be evaluated later with the
%   function nrbdeval2.
% 
% See also:
% 
%       nrbdeval
%
% Copyright (C) 2000 Mark Spink
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2010, 2011 Rafael Vazquez
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

if (~isstruct(nurbs))
  error('NURBS representation is not structure!');
end

if (~strcmp(nurbs.form,'B-NURBS'))
  error('Not a recognised NURBS representation');
end

% We raise the degree to avoid errors in the computation of the second derivative
if (nargout == 2)
  degelev  = max ([2 2] - (nurbs.order-1), 0);
  nurbs    = nrbdegelev (nurbs, degelev);
end

degree = nurbs.order - 1;

if (iscell(nurbs.knots))
  if (size(nurbs.knots,2) == 3)
  % NURBS structure represents a volume
    num1 = nurbs.number(1);
    num2 = nurbs.number(2);
    num3 = nurbs.number(3);

  % taking derivatives along the u direction
    dknots = nurbs.knots;
    dcoefs = permute (nurbs.coefs,[1 3 4 2]);
    dcoefs = reshape (dcoefs,4*num2*num3,num1);
    [dcoefs,dknots{1}] = bspderiv (degree(1),dcoefs,nurbs.knots{1});
    dcoefs = permute (reshape (dcoefs,[4 num2 num3 size(dcoefs,2)]),[1 4 2 3]);
    dnurbs{1} = nrbmak (dcoefs, dknots);

  % taking derivatives along the v direction
    dknots = nurbs.knots;
    dcoefs = permute (nurbs.coefs,[1 2 4 3]);
    dcoefs = reshape (dcoefs,4*num1*num3,num2);
    [dcoefs,dknots{2}] = bspderiv (degree(2),dcoefs,nurbs.knots{2});
    dcoefs = permute (reshape (dcoefs,[4 num1 num3 size(dcoefs,2)]),[1 2 4 3]);
    dnurbs{2} = nrbmak (dcoefs, dknots);

  % taking derivatives along the w direction
    dknots = nurbs.knots;
    dcoefs = reshape (nurbs.coefs,4*num1*num2,num3);
    [dcoefs,dknots{3}] = bspderiv (degree(3),dcoefs,nurbs.knots{3});
    dcoefs = reshape (dcoefs,[4 num1 num2 size(dcoefs,2)]);
    dnurbs{3} = nrbmak (dcoefs, dknots);

    if (nargout == 2)
      warning ('nrbderiv: the second derivative is not ready for volumes');
      dnurbs2 = [];
    end

  elseif (size(nurbs.knots,2) == 2)
% NURBS structure represents a surface

    num1 = nurbs.number(1);
    num2 = nurbs.number(2);

% taking first derivative along the u direction
    dknots = nurbs.knots;
    dcoefs = permute (nurbs.coefs,[1 3 2]);
    dcoefs = reshape (dcoefs,4*num2,num1);
    [dcoefs,dknots{1}] = bspderiv (degree(1),dcoefs,nurbs.knots{1});
    dcoefs = permute (reshape (dcoefs,[4 num2 size(dcoefs,2)]),[1 3 2]);
    dnurbs{1} = nrbmak (dcoefs, dknots);

    if (nargout == 2)
% taking second derivative along the u direction (duu)
      dknots2 = dknots;
      dcoefs2 = permute (dcoefs, [1 3 2]);
      dcoefs2 = reshape (dcoefs2, 4*num2, []);
      [dcoefs2, dknots2{1}] = bspderiv (degree(1)-1, dcoefs2, dknots{1});
      dcoefs2 = permute (reshape (dcoefs2, 4, num2, []), [1 3 2]);
      dnurbs2{1,1} = nrbmak (dcoefs2, dknots2); 

% taking second derivative along the v direction (duv and dvu)
      dknots2 = dknots;
      dcoefs2 = reshape (dcoefs, 4*(num1-1), num2);
      [dcoefs2, dknots2{2}] = bspderiv (degree(2), dcoefs2, dknots{2});
      dcoefs2 = reshape (dcoefs2, 4, num1-1, []);
      dnurbs2{1,2} = nrbmak (dcoefs2, dknots2);
      dnurbs2{2,1} = dnurbs2{1,2};
    end

% taking first derivative along the v direction
    dknots = nurbs.knots;
    dcoefs = reshape (nurbs.coefs,4*num1,num2);
    [dcoefs,dknots{2}] = bspderiv (degree(2),dcoefs,nurbs.knots{2});
    dcoefs = reshape (dcoefs,[4 num1 size(dcoefs,2)]);
    dnurbs{2} = nrbmak (dcoefs, dknots);

    if (nargout == 2)
% taking second derivative along the v direction (dvv)
      dknots2 = dknots;
      dcoefs2 = reshape (dcoefs, 4*num1, num2-1);
      [dcoefs2, dknots2{2}] = bspderiv (degree(2)-1, dcoefs2, dknots{2});
      dcoefs2 = reshape (dcoefs2, 4, num1, []);
      dnurbs2{2,2} = nrbmak (dcoefs2, dknots2);
    end

  end
else
  % NURBS structure represents a curve

  [dcoefs,dknots] = bspderiv (degree, nurbs.coefs, nurbs.knots);
  dnurbs = nrbmak (dcoefs, dknots);
  if (nargout == 2)
    [dcoefs2,dknots2] = bspderiv (degree-1, dcoefs, dknots);
    dnurbs2 = nrbmak (dcoefs2, dknots2);
  end

end

varargout{1} = dnurbs;
if (nargout == 2)
  varargout{2} = dnurbs2;
end

end

%!demo
%! crv = nrbtestcrv;
%! nrbplot(crv,48);
%! title('First derivatives along a test curve.');
%! 
%! tt = linspace(0.0,1.0,9);
%! 
%! dcrv = nrbderiv(crv);
%! 
%! [p1, dp] = nrbdeval(crv,dcrv,tt);
%! 
%! p2 = vecnorm(dp);
%! 
%! hold on;
%! plot(p1(1,:),p1(2,:),'ro');
%! h = quiver(p1(1,:),p1(2,:),p2(1,:),p2(2,:),0);
%! set(h,'Color','black');
%! hold off;

%!demo
%! srf = nrbtestsrf;
%! p = nrbeval(srf,{linspace(0.0,1.0,20) linspace(0.0,1.0,20)});
%! h = surf(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:)));
%! set(h,'FaceColor','blue','EdgeColor','blue');
%! title('First derivatives over a test surface.');
%!
%! npts = 5;
%! tt = linspace(0.0,1.0,npts);
%! dsrf = nrbderiv(srf);
%! 
%! [p1, dp] = nrbdeval(srf, dsrf, {tt, tt});
%! 
%! up2 = vecnorm(dp{1});
%! vp2 = vecnorm(dp{2});
%! 
%! hold on;
%! plot3(p1(1,:),p1(2,:),p1(3,:),'ro');
%! h1 = quiver3(p1(1,:),p1(2,:),p1(3,:),up2(1,:),up2(2,:),up2(3,:));
%! h2 = quiver3(p1(1,:),p1(2,:),p1(3,:),vp2(1,:),vp2(2,:),vp2(3,:));
%! set(h1,'Color','black');
%! set(h2,'Color','black');
%! 
%! hold off;
