function mnrb = nrbmodw (nrb, new_w, index)

%
% NRBMODW: Modify the weights of specific control points of any NURBS map.
%
% Calling Sequence:
% 
%   nrb = nrbmodw (nrb, new_w, index);
%   
%    INPUT:
%   
%      nrb   - NURBS map to be modified.
%      new_w - vector specifying the new weigths.
%      index - indices of the control points to be modified.
%   
%    OUTPUT:
%   
%      mnrb - the modified NURBS.
%   
% Copyright (C) 2015 Jacopo Corno
%
  
  mnrb = nrb;
  [ii, jj, kk] = ind2sub (nrb.number, index);
  for count = 1:numel (ii)
    mnrb.coefs(:,ii(count),jj(count),kk(count)) = ...
      [nrb.coefs(1:3,ii(count),jj(count),kk(count))./repmat(nrb.coefs(4,ii(count),jj(count),kk(count)),3,1).*new_w(count); new_w(count)];
  end
  
end

