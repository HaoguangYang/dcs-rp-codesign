%%%%%%%%%%% A CLI Tool for Distributed Controller %%%%%%%%%%%%
%%%%%%%%%%%%%% Resource / Performance Co-Design %%%%%%%%%%%%%%
%                                                            %
%   Closed Loop interface with sparcity tuning K generator   %
%   Two-pass First Fit Descend Consolidator                  %
%   Utilizes Ethernet for communication                      %
%   Author: Haoguang Yang                                    %
%   Date: 11/28/2017                                         %
%                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
close all;

max_num_of_processors = 17;
max_acceptable_performance_degrade_ratio = 0.05;
download_to_hardware = false;

% System Specification
load('humanInTheLoopTest.mat');

% hardware capacity profile
load('zynq_7000_spmv_capacity_profile.mat', 'wcet_table');
reserved_util = 0.30;

max_iter = 100;

unfulfilled_tasks = [];                      	% For initializing the loop

% Initialize
gamma = 0;
soln = H2sparse(A,Bd,B,Q,R,gamma);
K = soln.F;
n_tasks = size(B,2);
soln_dense_opt = soln;

iter = 1;
while (true)
    % Demultiplexor by rows
    blockified_tasks = Demultiplexor(K, n_tasks);

    % ET_Estimator, takes in blockified_tasks and HW capacity profile
    ET = ET_Estimator(blockified_tasks, wcet_table);

    % RMA Consolidator, Feed Back: tasks should be simplified a total factor
    % of shrink_factor (in terms of ET), to be fit into the computation nodes.
    [procs_occupied, tasks_on_proc, proc_util, proc_alloc_to_task, unfulfilled_tasks, shrink_factor] = ...
        RMA_Consolidator(max_num_of_processors, 1. - reserved_util, ET, freq);
    
    fprintf ('Used #Nodes: %i \t Tasks unalloc: %i \t Avg Node Util:%f\n', ...
        procs_occupied, length(unfulfilled_tasks), sum(proc_util(1:procs_occupied))/procs_occupied);

    if (~isempty(unfulfilled_tasks))
        % needs to further simplify
        gamma = gamma+1e-04;
        soln = H2sparse(A,Bd,B,Q,R,gamma);
        K = soln.F;
        % Check iteration boundary
        if(iter >= max_iter)
			error 'Error: Hardware constraints not met, max iterations reached!';
        end
        iter = iter + 1;
    else
        ratio = (soln.J - soln_dense_opt.J) / soln_dense_opt.J;
        if (ratio <= max_acceptable_performance_degrade_ratio)
            break;
        else
            fprintf ("Error: Hardware resource and Performance cannot be simultaneously satisfied:\n" + ...
                "\t With %i computing nodes the best achievable control performance is J=%f\n" + ...
                "\t This deviates from the dense optimal controller with J*=%f by %f.\n" + ...
                "What do you want to do:\n", ...
                procs_occupied, soln.J, soln_dense_opt.J, ratio);
            sel = input( ...
                "Y - increase the number of available nodes; \n" + ...
                "N - continue with the more degraded performance.\n", "s");
            if (sel == "Y" || sel == "y")
                max_num_of_processors = max_num_of_processors + 1
                iter = 1;
                gamma = 0;
                K = soln_dense_opt.F;
                soln = soln_dense_opt;
                continue;
            else
                break;
            end
        end
    end
end

disp('Verifying Solution...');
RMA_Verify(procs_occupied, tasks_on_proc, proc_util, proc_alloc_to_task, ET, 1. - reserved_util, freq);

if (~download_to_hardware)
    return;
end

% TRANSFER DATA TO HARDWARE
disp ('Initializing Hardware...');
%%%%%%%%%%%%%%%%% NETWORKING INTERFACE %%%%%%%%%%%%%%%%%%
% >>SEND
% 2*uint16_t            [n_tasks, NumPeriods]
% 2*uint32_t			bufferSize{send, recv};
% <<RECV
% 2*uint32_t			bufferSize{send, recv}; %Receive client and send host, two pass.
% >>SEND
% uint16_t              Task specific communication port number
%-----SWITCH TO SPECIFIC PORT-----
% >>SEND
% 2*uint16_t			size_A{nRow, nCol};
% nRow*uint32_t         A_ptr[1:size_A[0]]; %Non-zero element counter pointer, last element indicate total non-zero elements
% num_nonZero*uint32_t	A_ind;				%Index to non-zero elements in A_nonZero
% num_nonZero*double	A_nonZero;			%Non-zero elements in A.
% <<RECV NON-BLOCKED
% Computation Result    Loop.
%********************************************************
PubPort = 60001:(60000+max_num_of_processors);
TaskPort = 50001:(50000+n_tasks);
bufSize = uint32([0 0]);
for i = 1:max_num_of_processors
    PubComm(i) = tcpip('0.0.0.0',PubPort(i),...
                       'NetworkRole','server',...
                       'ByteOrder','littleEndian',...
                       'Timeout',10.0);
end
disp ('Waiting for Connection...');
for i = 1:max_num_of_processors
    fopen(PubComm(i));
    % Number of Public Port equals number of total nodes, each sending their corresponding task indexes.
    % Public port may only have buffer size of 512.
    network_stream_write(PubComm(i), uint16([length(tasks_on_proc{i}) 1]),512);    % Modify NumPeriods (1) here for future development.
    
    clientBufSize = network_stream_read(PubComm, 8, 512, 'uint32');
    % Task-specific ports have buffer size 65535.
    network_stream_write(PubComm, uint32([65535, 65535]),512);
    bufSize(1) = min(uint32(65535), clientBufSize(2));
    bufSize(2) = min(uint32(65535), clientBufSize(1));
end
% Task-specific port initialize
for i = 1:n_tasks
    % Task-specific ports have buffer size 65535.
    comm(i) = tcpip('0.0.0.0',TaskPort(i),...
                'NetworkRole','server',...
                'InputBufferSize',65535,...
                'OutputBufferSize',65535,...
                'ByteOrder','littleEndian',...
                'Timeout',30.0);
end

%Initialize each tasak
disp ('Launching tasks...');
for i = 1:n_tasks
    % Initiate on the corresponding node
    network_stream_write(PubComm(proc_alloc_to_task(i)), uint16(TaskPort(i)), bufSize(1));
    fopen(comm(i));
    BKEq = blockified_tasks{i, 1};
    shape = size(BKEq);
    size_BKEq = uint16(shape);
    network_stream_write(comm(i), size_BKEq, bufSize(1));
    
    % Write CSR Format, BEWARE C INDEX STARTS AT 0.
    [ptr, col, v] = CSR(BKEq);
    network_stream_write(comm(i), uint32(ptr), bufSize(1));
    col = col-1;
    network_stream_write(comm(i), uint32(col), bufSize(1));
    network_stream_write(comm(i), v, bufSize(1));
    % Write Dense matrix through socket
    % network_stream_write(comm(i), reshape(BKEq',1,[]), bufSize(1));
    fprintf ('Task %d initialized.\n',i);
end
% size of result from each task
nRow = zeros(n_tasks, 1);
for i = 1:n_tasks
    nRow(i) = size(blockified_tasks{i,1}, 1);
end

% Tasks run on node with release time = 0 and frequency = 25Hz, refer to
% QNX C Code for further modification.
for i = 1:max_num_of_processors
    % Start Trigger
    fwrite(PubComm(i), uint8(1), 'uint8');
end
disp('Collecting Computation Results...')
tmp = zeros(n_tasks, max(double(nRow(i))));
while(1)
    for i = 1:n_tasks
        %You may collect results here.
        tmp(i,1:nRow(i)) = fread(comm(i), double(nRow(i)), 'double');
    end
end
