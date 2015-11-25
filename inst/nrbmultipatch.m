function [interfaces, boundary] = nrbmultipatch (nurbs)

%
% NRBMULTIPATCH: construct the information for gluing conforming NURBS patches, using the same format as in GeoPDEs.
% 
% Calling Sequence:
% 
%   [interfaces, boundary] = nrbmultipatch (nurbs);
% 
% INPUT:
% 
%   nurbs   : an array of NURBS surfaces or volumes (not both), see nrbmak.
% 
% OUTPUT: 
% 
%   interfaces: array with the information for each interface, that is:
%      - number of the first patch (patch1), and the local side number (side1)
%      - number of the second patch (patch2), and the local side number (side2)
%      - flag (faces and volumes), ornt1, ornt2 (only volumes): information
%        on how the two patches match, see below.
%   boundary:   array with the boundary faces that do not belong to any interface
%      - nsides:  total number of sides on the boundary array (numel(boundary))
%      - patches: number of the patch to which the boundary belongs
%      - sides:   number of the local side on the patch
%
% The boundary faces are stored separately, that is, nsides=1 for each boundary.
%  To join several faces under the same condition, the user should do it by hand.
% 
%    Copyright (C) 2014, 2015 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


npatch = numel (nurbs);
if (~iscell (nurbs(1).knots))
  ndim = 1;
  face_corners = @(x) x.coefs;
  tol = 1e-15;
  compare_corners = @(nrb1, nrb2) max(abs(nrb1 - nrb2)) < tol;
elseif (size(nurbs(1).knots,2) == 2)
  ndim = 2;
  face_corners = @(x) x.coefs(:, [1 end]);
  compare_corners = @(nrb1, nrb2) compare_corners_univariate (nrb1, nrb2);
elseif (size(nurbs(1).knots,2) == 3)
  ndim = 3;
  face_corners = @(x) x.coefs(:, [1 end], [1 end]);
  compare_corners = @(nrb1, nrb2) compare_corners_bivariate (nrb1, nrb2);
end

non_set_faces = cell (npatch, 1);
for ii = 1:npatch
  if (~iscell (nurbs(ii).knots))
    if (ndim ~= 1)
      error ('All the patches must have the same dimension (at least for now)')
    end
  elseif (ndim ~= size(nurbs(ii).knots,2))
    error ('All the patches must have the same dimension (at least for now)')      
  end
  non_set_faces{ii} = 1:2*ndim;
end

num_interfaces = 0;

num_boundaries = 0;
boundary = struct ('nsides', 0, 'patches', [], 'faces', []);

for i1 = 1:npatch
  nrb_faces1 = nrbextract (nurbs(i1));
  for j1 = non_set_faces{i1}
    
% This is to fix a bug when two faces of the same patch form an interface
%  (for instance, in a ring or a torus)
    if (isempty (intersect (non_set_faces{i1}, j1))); continue; end

    nrb1 = nrb_faces1(j1);
    corners1 = face_corners (nrb1);

    non_set_faces{i1} = setdiff (non_set_faces{i1}, j1);
    flag = 0;

    i2 = i1 - 1;
    while (~flag && i2 < npatch)
      i2 = i2 + 1;
      nrb_faces2 = nrbextract (nurbs(i2));
      j2 = 0;
      while (~flag && j2 < numel (non_set_faces{i2}))
        j2 = j2 + 1;
        nrb2 = nrb_faces2(non_set_faces{i2}(j2));

        corners2 = face_corners (nrb2);
        if (ndim == 2 || ndim == 1)
          flag = compare_corners (corners1, corners2);
        elseif (ndim == 3)
          [flag, ornt1, ornt2] = compare_corners (corners1, corners2);
        end
        
        if (flag)
          if (numel(nrb1.coefs) ~= numel(nrb2.coefs))
            flag = 0;
            warning (['The corners of PATCH %d FACE %d, and PATCH %d FACE %d coincide, but the number ' ... 
                'of control points is different. No information is saved in this case'], i1, j1, i2, j2)
          end
        end
      end
    end

    if (flag)
      intrfc.patch1 = i1;
      intrfc.side1 = j1;
      intrfc.patch2 = i2;
      intrfc.side2 = non_set_faces{i2}(j2);
      if (ndim ==3)
        intrfc.flag = flag;
        intrfc.ornt1 = ornt1;
        intrfc.ornt2 = ornt2;
      elseif (ndim == 2)
        intrfc.ornt = flag;
      end

      non_set_faces{i2} = setdiff (non_set_faces{i2}, non_set_faces{i2}(j2));
      num_interfaces = num_interfaces + 1;
      interfaces(num_interfaces) = intrfc;
    else
      bndry.nsides = 1;
      bndry.patches = i1;
      bndry.faces = j1;
      num_boundaries = num_boundaries + 1;
      boundary(num_boundaries) = bndry;
    end
  end
end

if (num_interfaces == 0)
   interfaces = []; 
end
if (num_boundaries == 0)
  boundary = [];
end

end



function [flag, ornt1, ornt2] = compare_corners_bivariate (coefs1, coefs2)
  tol = 1e-13;

  coefs1 = reshape (coefs1, 4, []);
  coefs2 = reshape (coefs2, 4, []);
% Should use some sort of relative error
  if (max (max (abs (coefs1 - coefs2))) < tol)
    flag = 1; ornt1 = 1; ornt2 = 1;
  elseif (max (max (abs (coefs1 - coefs2(:,[1 3 2 4])))) < tol)
    flag = -1; ornt1 = 1; ornt2 = 1;
  elseif (max (max (abs (coefs1 - coefs2(:,[3 1 4 2])))) < tol)
    flag = -1; ornt1 = -1; ornt2 = 1;
  elseif (max (max (abs (coefs1 - coefs2(:,[2 1 4 3])))) < tol)
    flag = 1; ornt1 = -1; ornt2 = 1;
  elseif (max (max (abs (coefs1 - coefs2(:,[4 3 2 1])))) < tol)
    flag = 1; ornt1 = -1; ornt2 = -1;
  elseif (max (max (abs (coefs1 - coefs2(:,[4 2 3 1])))) < tol)
    flag = -1; ornt1 = -1; ornt2 = -1;
  elseif (max (max (abs (coefs1 - coefs2(:,[2 4 1 3])))) < tol)
    flag = -1; ornt1 = 1; ornt2 = -1;
  elseif (max (max (abs (coefs1 - coefs2(:,[3 4 1 2])))) < tol)
    flag = 1; ornt1 = 1; ornt2 = -1;
  else
    flag = 0; ornt1 = 0; ornt2 = 0;
  end
end

function flag = compare_corners_univariate (coefs1, coefs2)
  tol = 1e-13;

  coefs1 = reshape (coefs1, 4, []);
  coefs2 = reshape (coefs2, 4, []);
% Should use some sort of relative error
  if (max (max (abs (coefs1 - coefs2))) < tol)
    flag = 1;
  elseif (max (max (abs (coefs1 - coefs2(:,[end 1])))) < tol)
    flag = -1;
  else
    flag = 0;
  end
end
