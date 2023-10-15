% Best-Fit bin-packing algorithm.
% It allows the user to pass in a customizable bin capacity function, which
% further takes in current bin id, current assignment, and other user-provided
% arguments.
% Returns:
% bin_id    - ID of allocated bins for each item
% util      - bin utilization, this value is NOT normalized by bin_capacity_func.

function [bin_id, util] = bfd(v, bin_capacity_func, varargin)
    assert(all(v < 1., "all"));
    bin_id = zeros(size(v,1),1);
    util = zeros(size(v));
    if (isempty(v))
        return;
    elseif (nargin < 2)
        bin_capacity_func = @(x, bin_id, varargin) 1.;
    end
    [v_sort, v_rank] = custom_sort(v, @(x) x(:,1), "descend");
    occupied_bin = 1;
    util(1,:) = v_sort(1, :);
    bin_id(v_rank(1)) = 1;
    bin_rank = (1);
    for i = 2 : size(v, 1)
        for j = 1 : occupied_bin
            util_tmp = util(bin_rank(j), :) + v_sort(i, :);
            bin_id(v_rank(i)) = j;
            if (any(util_tmp > bin_capacity_func(j, bin_id, varargin), "all"))
                if (j == occupied_bin)
                    occupied_bin = occupied_bin + 1;
                    util(occupied_bin, :) = v_sort(i, :);
                    bin_id(v_rank(i)) = occupied_bin;
                    [~, bin_rank] = custom_sort(util(1:occupied_bin, :), @(x) x(:,1), "descend");
                    break;
                else
                    bin_id(v_rank(i)) = 0;
                    continue;
                end
            else
                util(bin_rank(j), :) = util_tmp;
                [~, bin_rank] = custom_sort(util(1:occupied_bin, :), @(x) x(:,1), "descend");
                break;
            end
        end
    end
    util = util(1:occupied_bin, :);
end
