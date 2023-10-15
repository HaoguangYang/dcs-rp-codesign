% A test case of the task consolidator using contrived task sets

clear all;
clc;
max_num_of_nodes = 100;
used_num_of_nodes = 0;
node_util_lim = 1;
loop = 200;

% task_utilization efficiency testing
% mean utilization of the arbitrated set of tasks.
mu = 0.15;
% mu_array = 0.1:0.01:0.9;
% for j = 1:size(mu,2)
%     mu = mu_array(j);
% stdev of the utilization
sigma = 0.075;

% Complexity testing
num_tasks = 100;
% num_tasks_array = [10:10:100, 125:25:500, 550:50:1000, 1100:100:1800];
% for j = 1:size(num_tasks_array, 2)
%     num_tasks = num_tasks_array(j);

    for i = 1:loop
        freq_multiplier = ceil(5 * rand(num_tasks,1));

        task_util = [];
        while (num_tasks-length(task_util) > 0)
            task_util = [task_util; random('Normal', mu, sigma, num_tasks-length(task_util), 1)];
            task_util = task_util(task_util>0);
            task_util = task_util(task_util<1);
        end
        ET = task_util./freq_multiplier;
        
        tic();
        
        [used_num_of_nodes, NumAssignedtoNode, node_util, IndicatorAssignedtoNode, ~, ~] = ...
           RMA_Consolidator(max_num_of_nodes, node_util_lim, ET, freq_multiplier);

        %[used_num_of_nodes, NumAssignedtoNode, node_util, IndicatorAssignedtoNode, ~, ~, ~] = ...
        %    OnePassFFD_Consolidator(max_num_of_nodes, node_util_lim, ET, freq_multiplier, false);

        %[used_num_of_nodes, NumAssignedtoNode, node_util, IndicatorAssignedtoNode, ~, ~, ~] = ...
        %    RMASequential_Consolidator(max_num_of_nodes, node_util_lim, ET, freq_multiplier, false);
        
        timing(i) = toc();
        
        disp ('Number of Nodes Used:');
        disp (used_num_of_nodes);

        util(i) = sum(node_util(1:used_num_of_nodes))/used_num_of_nodes;
        disp ('Average Node Utilization:');
        disp (util(i));

        nnode(i) = used_num_of_nodes;
        assnability(i) = 1-sum(IndicatorAssignedtoNode==0)/num_tasks;
        
        disp('Verifying Solution...');
        RMA_Verify(used_num_of_nodes, NumAssignedtoNode, node_util, IndicatorAssignedtoNode, ET, node_util_lim, freq_multiplier);
    end

%    avg_time(j) = sum(timing)/loop;
%    avg_util(j) = sum(util)/loop;
%    avg_assn(j) = sum(assnability)/loop;
%    wc_assn(j) = min(assnability);

%end
