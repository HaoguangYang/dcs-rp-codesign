% calling sort using a user-defined function fcn.

function [sorted, rank] = custom_sort(elems, fcn, direction)
    elem_cost = fcn(elems);
    [~, rank] = sort(elem_cost, direction);
    sorted = elems(rank, :);
end
