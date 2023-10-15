% In a bin-packing setting, prune bins with utilization lower than the lower_bound.
% Bin X+1 following the pruned bin X will be re-assigned a new ID of X.
% original items in bin X will be labled as unassigned (bin_id = 0).

function [bin_id, util] = prune_loose_bins(bin_id_in, util_in, lower_bound)
    bin_id = bin_id_in;
    util = util_in;
    curr_bin = 1;
    for i = 1 : size(util, 1) - 1
        if (util(curr_bin) < lower_bound)
            % migrate items from bin #i+1 to bin #curr_bin
            util(curr_bin) = util(i + 1);
            bin_id(bin_id == i) = 0;
        else
            % i <- i + 1 in for...
            bin_id(bin_id == i) = curr_bin;
            curr_bin = curr_bin + 1;
        end
    end
    % last loop
    if (util(curr_bin) < lower_bound)
        bin_id(bin_id == size(util, 1)) = 0;
    else
        % i <- i + 1 in for...
        bin_id(bin_id == size(util, 1)) = curr_bin;
        curr_bin = curr_bin + 1;
    end
    util = util(1 : curr_bin - 1);
end
