function [pnt,jac] = nrbdeval(nurbs, dnurbs, tt)

% NRBDEVAL: Evaluation of the derivative NURBS curve, surface or volume.
%
%     [pnt, jac] = nrbdeval(crv, dcrv, tt)
%     [pnt, jac] = nrbdeval(srf, dsrf, {tu tv})
%     [pnt, jac] = nrbdeval(vol, dvol, {tu tv tw})
%
% INPUTS:
%
%   crv    - original NURBS curve.
%
%   srf    - original NUBRS surface
%
%   vol    - original NUBRS volume
%
%   dcrv   - NURBS derivative representation of crv
%
%   dsrf   - NURBS derivative representation of surface
%
%   dvol   - NURBS derivative representation of volume
%
%   tt     - parametric evaluation points
%            If the nurbs is a surface or a volume then tt is a cell
%            {tu, tv} or {tu, tv, tw} are the parametric coordinates
%
% OUTPUT:
%
%   pnt  - evaluated points.
%   jac  - evaluated first derivatives (Jacobian).
%
% Examples:
% 
%   // Determine the first derivatives a NURBS curve at 9 points for 0.0 to
%   // 1.0
%   tt = linspace(0.0, 1.0, 9);
%   dcrv = nrbderiv(crv);
%   [pnts,jac] = nrbdeval(crv, dcrv, tt);
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

if (~isstruct(nurbs))
  error('NURBS representation is not structure!');
end

if (~strcmp(nurbs.form,'B-NURBS'))
  error('Not a recognised NURBS representation');
end

[cp,cw] = nrbeval (nurbs, tt);

if (iscell(nurbs.knots))
  if (size(nurbs.knots,2) == 3)
  % NURBS structure represents a volume
    temp = cw(ones(3,1),:,:,:);
    pnt = cp./temp;
  
    [cup,cuw] = nrbeval (dnurbs{1}, tt);
    tempu = cuw(ones(3,1),:,:,:);
    jac{1} = (cup-tempu.*pnt)./temp;
  
    [cvp,cvw] = nrbeval (dnurbs{2}, tt);
    tempv = cvw(ones(3,1),:,:,:);
    jac{2} = (cvp-tempv.*pnt)./temp;

    [cwp,cww] = nrbeval (dnurbs{3}, tt);
    tempw = cww(ones(3,1),:,:,:);
    jac{3} = (cwp-tempw.*pnt)./temp;

  elseif (size(nurbs.knots,2) == 2)
  % NURBS structure represents a surface
    temp = cw(ones(3,1),:,:);
    pnt = cp./temp;
  
    [cup,cuw] = nrbeval (dnurbs{1}, tt);
    tempu = cuw(ones(3,1),:,:);
    jac{1} = (cup-tempu.*pnt)./temp;
  
    [cvp,cvw] = nrbeval (dnurbs{2}, tt);
    tempv = cvw(ones(3,1),:,:);
    jac{2} = (cvp-tempv.*pnt)./temp;

  end
else

  % NURBS is a curve
  temp = cw(ones(3,1),:);
  pnt = cp./temp;
  
  % first derivative
  [cup,cuw] = nrbeval (dnurbs,tt);
  temp1 = cuw(ones(3,1),:);
  jac = (cup-temp1.*pnt)./temp;

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

%!test
%! knots{1} = [0 0 0 1 1 1];
%! knots{2} = [0 0 0 .5 1 1 1];
%! knots{3} = [0 0 0 0 1 1 1 1];
%! cx = [0 0.5 1]; nx = length(cx);
%! cy = [0 0.25 0.75 1]; ny = length(cy);
%! cz = [0 1/3 2/3 1]; nz = length(cz);
%! coefs(1,:,:,:) = repmat(reshape(cx,nx,1,1),[1 ny nz]);
%! coefs(2,:,:,:) = repmat(reshape(cy,1,ny,1),[nx 1 nz]);
%! coefs(3,:,:,:) = repmat(reshape(cz,1,1,nz),[nx ny 1]);
%! coefs(4,:,:,:) = 1;
%! nurbs = nrbmak(coefs, knots);
%! x = rand(5,1); y = rand(5,1); z = rand(5,1);
%! tt = [x y z]';
%! ders = nrbderiv(nurbs);
%! [points,jac] = nrbdeval(nurbs,ders,tt);
%! assert(points,tt,1e-10)
%! assert(jac{1}(1,:,:),ones(size(jac{1}(1,:,:))),1e-12)
%! assert(jac{2}(2,:,:),ones(size(jac{2}(2,:,:))),1e-12)
%! assert(jac{3}(3,:,:),ones(size(jac{3}(3,:,:))),1e-12)
%! 
%!test
%! knots{1} = [0 0 0 1 1 1];
%! knots{2} = [0 0 0 0 1 1 1 1];
%! knots{3} = [0 0 0 1 1 1];
%! cx = [0 0 1]; nx = length(cx);
%! cy = [0 0 0 1]; ny = length(cy);
%! cz = [0 0.5 1]; nz = length(cz);
%! coefs(1,:,:,:) = repmat(reshape(cx,nx,1,1),[1 ny nz]);
%! coefs(2,:,:,:) = repmat(reshape(cy,1,ny,1),[nx 1 nz]);
%! coefs(3,:,:,:) = repmat(reshape(cz,1,1,nz),[nx ny 1]);
%! coefs(4,:,:,:) = 1;
%! coefs = coefs([2 1 3 4],:,:,:);
%! nurbs = nrbmak(coefs, knots);
%! x = rand(5,1); y = rand(5,1); z = rand(5,1);
%! tt = [x y z]';
%! dnurbs = nrbderiv(nurbs);
%! [points, jac] = nrbdeval(nurbs,dnurbs,tt);
%! assert(points,[y.^3 x.^2 z]',1e-10);
%! assert(jac{2}(1,:,:),3*y'.^2,1e-12)
%! assert(jac{1}(2,:,:),2*x',1e-12)
%! assert(jac{3}(3,:,:),ones(size(z')),1e-12)