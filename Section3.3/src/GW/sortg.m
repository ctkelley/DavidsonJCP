function sgidx = sortg(gmask)
%
% sgidx = sortg(gmask)
% 
% Sort the reciprocal space indices according to their distance
% to the original in ascending order. Indices that have the same 
% distances to the origin are arranged in lexicographical order
%
[sgkk, sgidx] = sort(gmask.gkk);
