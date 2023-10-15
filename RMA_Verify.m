% A function to verify the result of Multi-Processor RMA Consolidator.

function RMA_Verify(occupied_procs, tasks_on_each_proc, util_of_proc, ...
                    proc_alloc_to_task, ET, max_et_avail_per_sec, freq)
    [n_tasks, ~] = size(ET);
    task_util = ET .* freq / max_et_avail_per_sec;
	
	is_allocation_valid = zeros(n_tasks, 1);
	unique_periods_on_proc = zeros(occupied_procs, 1);
	correct = 1;
	for i = 1 : occupied_procs
		j = 1;
        task_on_proc = tasks_on_each_proc{i};
        while (j <= length(task_on_proc) && proc_alloc_to_task(task_on_proc(j)) == i)
			is_allocation_valid(task_on_proc(j)) = 1;
			j = j + 1;
        end
        unique_periods_on_proc(i) = length(unique(freq(task_on_proc)));
		util_lim_rm_sched = unique_periods_on_proc(i) * (2 ^ (1. / unique_periods_on_proc(i)) - 1.);
        if (abs(sum(task_util(task_on_proc), 1) - util_of_proc(i)) > 1e-6)
            disp ('Node Util Calculation is Incorrect!');
			disp (i);
            correct = 0;
        end
		if util_of_proc(i) > util_lim_rm_sched
			disp ('Node Util Exceed Limit!');
			disp (i);
			correct = 0;
		end
	end
	if (~all(is_allocation_valid))
		disp ('Tasks not corectly allocated!');
		correct = 0;
	end
	if (correct)
        disp ('Solution is Correct!');
	else
		disp ('Solution is Incorrect!');
	end
end
