% WCET estimator
% This component requires a hardware capacity profile obtained by profiling
% a set of operation stencils on the target hardware.
% example file of the hardware capacity profile can be generated with:
%
% % number of stencils (rows) * numer of operations (cols)
% wcet_table = cell(2,2);
% % operation stencils: scale function w.r.t. computing input x
% wcet_table{1, 1} = @(x) sum(x ~= 0, 'all');
% wcet_table{2, 1} = @(x) 1.;
% % unit stencil WCETs
% wcet_table{1, 2} = 4.133e-8;
% wcet_table{2, 2} = 2.658e-5;
% save("zynq_7000_spmv_capacity_profile.mat", "wcet_table");

function ET = ET_Estimator(computing_tasks, wcet_table)
    n_tasks = size(computing_tasks, 1);
    ET = zeros(n_tasks, 1);
    for i = 1 : n_tasks
       for j = 1 : size(wcet_table, 1)
           ET(i,1) = ET(i,1) + wcet_table{j,1}(computing_tasks{i, 1}) * wcet_table{j,2};
       end
    end
end
