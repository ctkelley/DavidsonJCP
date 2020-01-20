function [idx] = cell_union(C)
idx = C{1};
for k = 2:length(C)
    idx = union(idx,C{k});
end