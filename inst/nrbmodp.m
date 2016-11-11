function mnrb = nrbmodp (nrb, move, index)

%
% NRBMODP: Modify the coordinates of specific control points of any NURBS
% map. The weight is not changed.
%
% Calling Sequence:
% 
%   nrb = nrbmodp (nrb, move, index);
%   
%    INPUT:
%   
%      nrb   - NURBS map to be modified.
%      move  - vector specifying the displacement of all the ctrl points.
%      index - indeces of the control points to be modified.
%   
%    OUTPUT:
%   
%      mnrb - the modified NURBS.
%   
% Copyright (C) 2015 Jacopo Corno
%
  
  move = reshape (move, 3, 1);

  mnrb = nrb;
  [ii, jj, kk] = ind2sub (nrb.number, index);
  for count = 1:numel (ii)
    mnrb.coefs(1:3,ii(count),jj(count),kk(count)) = nrb.coefs(1:3,ii(count),jj(count),kk(count)) + ...
      move * nrb.coefs(4,ii(count),jj(count),kk(count));
  end
  
end

