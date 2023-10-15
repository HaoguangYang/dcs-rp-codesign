% calling sort using a user-defined function fcn.

function util_lim = rm_sched_util_limit(id, all_allocs, task_periods)
    n_unique_periods = length(unique(task_periods(all_allocs == id)));
    util_lim = n_unique_periods * (2 ^ (1./n_unique_periods) - 1.);
end
