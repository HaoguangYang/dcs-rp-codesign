% Multi-Processor RMA Consolidator based on 2-Pass FFD Algorithm
% Author: 			Haoguang Yang
% Last Modified:	Jan-19-2018
% Algorithm explained on Section VII of the GE Report.

function [proc_occupied, tasks_on_proc, proc_util, proc_allocated, ...
            tasks_to_shrink, shrink_factor] = ...
            RMA_Consolidator(proc_lim, max_et_avail_per_sec, ET, freq)
    
    quiet = 1;
    % Multi-Processor RMA Consolidator
    % Allocation
    [~, order] = sort(freq);
    [n_tasks, ~] = size(ET);
    et_per_second = ET .* freq;
    task_util = et_per_second / max_et_avail_per_sec;
    total_exec_time = sum(et_per_second);

    shrink_factor = 0.;
    tasks_to_shrink(1,:) = find(task_util > 1.);
    % exists a single task that overruns capacity and needs further simplification.
    if (~isempty(tasks_to_shrink))
        shrink_factor = max(task_util(tasks_to_shrink),1) - 1.;
        shrink_factor(shrink_factor < 0.) = 0.;
        error 'Task Exceeded Node Computational Limit!';
    end
	
	period_group = zeros(n_tasks,1);    % How many tasks in corresponding periods 
	n_unique_periods = 1;               % Number of Different Periods in the task pool
    period_group(1) = 1;
    for i = 2 : n_tasks
        if (freq(order(i)) == freq(order(i-1)))
			period_group(n_unique_periods) = period_group(n_unique_periods) + 1;
		else
			n_unique_periods = n_unique_periods + 1;
			period_group(n_unique_periods) = 1;
        end
    end
	
    util_lim_one_proc = n_unique_periods * (2 ^ (1 / n_unique_periods) - 1);
    
    % obtain a rough estimate of how many nodes needed
    exec_time_per_node_est = max_et_avail_per_sec * util_lim_one_proc;
    % [quot, rem] = quorem(sym(total_exec_time), sym(max_et_avail_per_sec));
	quot = floor(total_exec_time / exec_time_per_node_est);
	rem = mod(total_exec_time, exec_time_per_node_est);
    if(rem > 0)
        proc_n_est = uint32(quot) + 1;
    else
        proc_n_est = uint32(quot);
    end
    
    if (~quiet)
        disp('Expected Number of Nodes Used:');
        disp(proc_n_est);
    end
    
    % Utilization of each processor
	proc_util = zeros(n_tasks, 1);
    % On which nodes is the task assigned
    proc_allocated = zeros(n_tasks, 1);
    % Number of nodes utilized
    proc_occupied = 0;
	
	for i = 1 : n_unique_periods    % Group tasks with same periods.
		period_group_base = sum(period_group(1 : i - 1)) + 1;
        if (period_group(i)==2)     % 2 grouped tasks would fit in 1 node.
			util_tmp = task_util(order(period_group_base)) + task_util(order(period_group_base + 1));
            if (util_tmp < 1 && util_tmp > 2 * (sqrt(2) - 1) && proc_occupied < proc_lim)
				% The two tasks fits in one node and sufficiently occupied it. Assign tasks to node
                proc_occupied = proc_occupied + 1;
				proc_allocated(order(period_group_base)) = proc_occupied;
				proc_allocated(order(period_group_base + 1)) = proc_occupied;
				% Update node utilitization
				proc_util(proc_occupied) = util_tmp;
            end
        elseif (period_group(i)>2)	% Use FFD as packing algorithm for task groups more than 3.
            task_id_in_group = order(period_group_base : period_group_base + period_group(i) - 1);
            %% FFD bin-packing method
            [proc_id_tmp, proc_util_tmp] = ffd(task_util(task_id_in_group));

            %% prune loose bins and re-arrange results
            [proc_id_tmp, proc_util_tmp] = prune_loose_bins(proc_id_tmp, proc_util_tmp, 2 * (sqrt(2) - 1));
		    proc_util(proc_occupied + 1 : proc_occupied + size(proc_util_tmp,1)) = proc_util_tmp;
            proc_id_tmp(proc_id_tmp > 0) = proc_id_tmp(proc_id_tmp > 0) + proc_occupied;
            proc_allocated(task_id_in_group) = proc_id_tmp;
            proc_occupied = proc_occupied + size(proc_util_tmp,1);
        end
	end
	
	% for remaining tasks, use general FFD algorithms
	remaining_tasks = find(proc_allocated == 0);

    %% FFD bin-packing method
    [proc_id_tmp, proc_util_tmp] = ffd(task_util(remaining_tasks), ...
        @(cur_proc, alloc_plan, periods_cell) rm_sched_util_limit(cur_proc, alloc_plan, periods_cell{1}), ...
        freq(remaining_tasks));
    
    %% combine results
    proc_allocated(remaining_tasks) = proc_id_tmp + proc_occupied;
    if (proc_occupied + size(proc_util_tmp,1) > proc_lim)
        fulfilled = proc_lim - proc_occupied;
        proc_util(proc_occupied + 1 : proc_lim) = proc_util_tmp(1 : fulfilled);
        tasks_to_shrink = find(proc_allocated > proc_lim);
        proc_allocated(tasks_to_shrink) = 0;
        shrink_factor = sum(task_util(tasks_to_shrink, :), 1)./sum(task_util, 1);
    else
        proc_util(proc_occupied + 1 : proc_occupied + size(proc_util_tmp,1)) = proc_util_tmp;
        proc_occupied = proc_occupied + size(proc_util_tmp,1);
    end
    
    tasks_on_proc = cell(proc_occupied, 1);
    for i = 1 : n_tasks
        tasks_on_proc{proc_allocated(i)}(end + 1) = i;
    end
    
    proc_util = proc_util(1:proc_occupied);
    % Scheduling
    % Scheduling of tasks is not included in this scope.
end
