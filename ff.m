% First-Fit bin-packing algorithm.
% It allows the user to pass in a customizable bin capacity function, which
% further takes in current bin id, current assignment, and other user-provided
% arguments.
% Returns:
% bin_id    - ID of allocated bins for each item
% util      - bin utilization, this value is NOT normalized by bin_capacity_func.

function [bin_id, util] = ff(v, bin_capacity_func, varargin)
    assert(all(v < 1., "all"));
    bin_id = zeros(size(v,1),1);
    util = zeros(size(v));
    if (isempty(v))
        return;
    elseif (nargin < 2)
        bin_capacity_func = @(x, bin_id, varargin) 1.;
    end
    occupied_bin = 1;
    util(1,:) = v(1, :);
    bin_id(1) = 1;
    for i = 2:size(v, 1)
        for j = 1:occupied_bin
            util_tmp = util(j, :) + v(i, :);
            bin_id(i) = j;
            if (any(util_tmp > bin_capacity_func(j, bin_id, varargin), "all"))
                if (j == occupied_bin)
                    occupied_bin = occupied_bin + 1;
                    util(occupied_bin, :) = v(i, :);
                    bin_id(i) = occupied_bin;
                    break;
                else
                    bin_id(i) = 0;
                    continue;
                end
            else
                util(j, :) = util_tmp;
                break;
            end
        end
    end
    util = util(1:occupied_bin, :);
end
