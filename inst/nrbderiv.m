function dnurbs = nrbderiv(nurbs)
% 
% NRBDERIV: Construct the first derivative representation of a
%           NURBS curve, surface or volume.
% 
% Calling Sequence:
% 
%   ders = nrbderiv(nrb);
% 
% INPUT:
% 
%   nrb		: NURBS data structure, see nrbmak.
%
% OUTPUT:
% 
%   ders	: A data structure that represents the first
% 		  derivatives of a NURBS curve, surface or volume.
% 
% Description:
% 
%   The derivatives of a B-Spline are themselves a B-Spline of lower degree,
%   giving an efficient means of evaluating multiple derivatives. However,
%   although the same approach can be applied to NURBS, the situation for
%   NURBS is a more complex. I have at present restricted the derivatives
%   to just the first. I don't claim that this implentation is
%   the best approach, but it will have to do for now. The function returns
%   a data struture that for a NURBS curve contains the first derivatives of
%   the B-Spline representation. Remember that a NURBS curve is represent by
%   a univariate B-Spline using the homogeneous coordinates.
%   The derivative data structure can be evaluated later with the function
%   nrbdeval.
% 
% Examples:
% 
%   See the function nrbdeval for an example.
% 
% See also:
% 
%       nrbdeval
%
% Copyright (C) 2000 Mark Spink
% Copyright (C) 2010 Rafael Vazquez
% Copyright (C) 2010 Carlo de Falco
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

if ~isstruct(nurbs)
  error('NURBS representation is not structure!');
end

if ~strcmp(nurbs.form,'B-NURBS')
  error('Not a recognised NURBS representation');
end

degree = nurbs.order - 1;

if iscell(nurbs.knots)
  if size(nurbs.knots,2) == 3
  % NURBS structure represents a volume
    num1 = nurbs.number(1);
    num2 = nurbs.number(2);
    num3 = nurbs.number(3);

  % taking derivatives along the u direction
    dknots = nurbs.knots;
    dcoefs = permute(nurbs.coefs,[1 3 4 2]);
    dcoefs = reshape(dcoefs,4*num2*num3,num1);
    [dcoefs,dknots{1}] = bspderiv(degree(1),dcoefs,nurbs.knots{1});
    dcoefs = permute(reshape(dcoefs,[4 num2 num3 size(dcoefs,2)]),[1 4 2 3]);
    dnurbs{1} = nrbmak(dcoefs, dknots);

  % taking derivatives along the v direction
    dknots = nurbs.knots;
    dcoefs = permute(nurbs.coefs,[1 2 4 3]);
    dcoefs = reshape(dcoefs,4*num1*num3,num2);
    [dcoefs,dknots{2}] = bspderiv(degree(2),dcoefs,nurbs.knots{2});
    dcoefs = permute(reshape(dcoefs,[4 num1 num3 size(dcoefs,2)]),[1 2 4 3]);
    dnurbs{2} = nrbmak(dcoefs, dknots);

  % taking derivatives along the w direction
    dknots = nurbs.knots;
    dcoefs = reshape(nurbs.coefs,4*num1*num2,num3);
    [dcoefs,dknots{3}] = bspderiv(degree(3),dcoefs,nurbs.knots{3});
    dcoefs = reshape(dcoefs,[4 num1 num2 size(dcoefs,2)]);
    dnurbs{3} = nrbmak(dcoefs, dknots);

  elseif size(nurbs.knots,2) == 2
  % NURBS structure represents a surface

    num1 = nurbs.number(1);
    num2 = nurbs.number(2);

  % taking derivatives along the u direction
    dknots = nurbs.knots;
    dcoefs = permute(nurbs.coefs,[1 3 2]);
    dcoefs = reshape(dcoefs,4*num2,num1);
    [dcoefs,dknots{1}] = bspderiv(degree(1),dcoefs,nurbs.knots{1});
    dcoefs = permute(reshape(dcoefs,[4 num2 size(dcoefs,2)]),[1 3 2]);
    dnurbs{1} = nrbmak(dcoefs, dknots);

  % taking derivatives along the v direction
    dknots = nurbs.knots;
    dcoefs = reshape(nurbs.coefs,4*num1,num2);
    [dcoefs,dknots{2}] = bspderiv(degree(2),dcoefs,nurbs.knots{2});
    dcoefs = reshape(dcoefs,[4 num1 size(dcoefs,2)]);
    dnurbs{2} = nrbmak(dcoefs, dknots);
  end
else
  % NURBS structure represents a curve

  [dcoefs,dknots] = bspderiv(degree,nurbs.coefs,nurbs.knots);
  dnurbs = nrbmak(dcoefs, dknots);

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
