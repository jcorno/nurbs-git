function nrbplot(nurbs,subd,p1,v1)
% 
% NRBPLOT: Plot a NURBS curve or surface.
% 
% Calling Sequence:
% 
%   nrbplot(nrb,subd)
%   nrbplot(nrb,subd,p,v)
% 
% INPUT:
% 
%   nrb		: NURBS curve or surface, see nrbmak.
% 
%   npnts	: Number of evaluation points, for a surface a row vector
% 		with two elements for the number of points along the U and
% 		V directions respectively.
% 
%   [p,v]       : property/value options
%
%               Valid property/value pairs include:
%
%               Property        Value/{Default}
%               -----------------------------------
%               light           {off} | true  
%               colormap        {'copper'}
%
% Example:
%
%   Plot the test surface with 20 points along the U direction
%   and 30 along the V direction
%
%   plot(nrbtestsrf, [20 30])
%
%    Copyright (C) 2000 Mark Spink
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

nargs = nargin;
if nargs < 2
  error('Need a NURBS to plot and the number of subdivisions!');
elseif rem(nargs+2,2)
  error('Param value pairs expected')
end

% Default values
light='off';
cmap='copper';

% Recover Param/Value pairs from argument list
for i=3:2:nargs
  Param = eval(['p' int2str((i-3)/2 +1)]);
  Value = eval(['v' int2str((i-3)/2 +1)]);
  if ~ischar(Param)
    error('Parameter must be a string')
  elseif size(Param,1)~=1
    error('Parameter must be a non-empty single row string.')
  end
  switch lower(Param)
  case 'light'
    light = lower(Value);
    if ~ischar(light)
      error('light must be a string.')
    elseif ~(strcmp(light,'off') | strcmp(light,'on'))
      error('light must be off | on')
    end
  case 'colormap'
    if ischar(Value)
      cmap = lower(Value);
    elseif size(Value,2) ~= 3
      error('colormap must be a string or have exactly three columns.')
    else
      cmap=Value;
    end
  otherwise
    error(['Unknown parameter: ' Param])
  end
end

colormap('default');

% convert the number of subdivisions in number of points
subd = subd+1;

% plot the curve or surface
if iscell(nurbs.knots)
 if size(nurbs.knots,2) == 2 % plot a NURBS surface
  p = nrbeval(nurbs,{linspace(0.0,1.0,subd(1)) linspace(0.0,1.0,subd(2))});
  if strcmp(light,'on') 
    % light surface
    surfl(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:)));
    shading interp;
    colormap(cmap);
  else 
    surf(squeeze(p(1,:,:)),squeeze(p(2,:,:)),squeeze(p(3,:,:)));
    shading faceted;
  end
 else
  error('The function nrbplot is not yet ready for volumes')
 end
else
  % plot a NURBS curve
  p = nrbeval(nurbs,linspace(0.0,1.0,subd));

  if any(nurbs.coefs(3,:))
    % 3D curve
    plot3(p(1,:),p(2,:),p(3,:)); 
    grid on;
  else
    % 2D curve
    plot(p(1,:),p(2,:));
  end
end
axis equal;

end

% plot the control surface
% hold on;
% mesh(squeeze(pnts(1,:,:)),squeeze(pnts(2,:,:)),squeeze(pnts(3,:,:)));
% hold off;

%!demo
%! crv = nrbtestcrv;
%! nrbplot(crv,100)
%! title('Test curve')
%! hold off

%!demo
%! coefs = [0.0 7.5 15.0 25.0 35.0 30.0 27.5 30.0;
%!          0.0 2.5  0.0 -5.0  5.0 15.0 22.5 30.0];
%! knots = [0.0 0.0 0.0 1/6 1/3 1/2 2/3 5/6 1.0 1.0 1.0];
%!
%! geom = [
%! nrbmak(coefs,knots)
%! nrbline([30.0 30.0],[20.0 30.0])
%! nrbline([20.0 30.0],[20.0 20.0])
%! nrbcirc(10.0,[10.0 20.0],1.5*pi,0.0)
%! nrbline([10.0 10.0],[0.0 10.0])
%! nrbline([0.0 10.0],[0.0 0.0])
%! nrbcirc(5.0,[22.5 7.5])
%! ];
%!
%! ng = length(geom);
%! for i = 1:ng
%!   nrbplot(geom(i),500);
%!   hold on;
%! end
%! hold off;
%! axis equal;
%! title('2D Geometry formed by a series of NURBS curves');

%!demo
%! sphere = nrbrevolve(nrbcirc(1,[],0.0,pi),[0.0 0.0 0.0],[1.0 0.0 0.0]);
%! nrbplot(sphere,[40 40],'light','on');
%! title('Ball and torus - surface construction by revolution');
%! hold on;
%! torus = nrbrevolve(nrbcirc(0.2,[0.9 1.0]),[0.0 0.0 0.0],[1.0 0.0 0.0]);
%! nrbplot(torus,[40 40],'light','on');
%! hold off
