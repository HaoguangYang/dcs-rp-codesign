% GE Control Sparsity Tuner and Demo Run
% Author: 			Haoguang Yang
% Last Modified:	Jan-19-2018
% Based on Dohyeung Kim's original GUI design.

function varargout = DistributedControllerGenerator(varargin)
% DISTRIBUTEDCONTROLLERGENERATOR MATLAB code for DistributedControllerGenerator.fig
%      DISTRIBUTEDCONTROLLERGENERATOR, by itself, creates a new DISTRIBUTEDCONTROLLERGENERATOR or raises the existing
%      singleton*.
%
%      H = DISTRIBUTEDCONTROLLERGENERATOR returns the handle to a new DISTRIBUTEDCONTROLLERGENERATOR or the handle to
%      the existing singleton*.
%
%      DISTRIBUTEDCONTROLLERGENERATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISTRIBUTEDCONTROLLERGENERATOR.M with the given input arguments.
%
%      DISTRIBUTEDCONTROLLERGENERATOR('Property','Value',...) creates a new DISTRIBUTEDCONTROLLERGENERATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DistributedControllerGenerator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DistributedControllerGenerator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DistributedControllerGenerator

% Last Modified by GUIDE v2.5 21-May-2020 11:37:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DistributedControllerGenerator_OpeningFcn, ...
                   'gui_OutputFcn',  @DistributedControllerGenerator_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DistributedControllerGenerator is made visible.
function DistributedControllerGenerator_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = DistributedControllerGenerator_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function NOACEdit_Callback(hObject, eventdata, handles)
data = guidata(hObject);
data.proc_lim = str2num(get(handles.NOACEdit, 'String'));
guidata(hObject, data);

% --- Executes during object creation, after setting all properties.
function NOACEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CDCEdit_Callback(hObject, eventdata, handles)
data = guidata(hObject);
data.ContDutyCycle = str2num(get(handles.CDCEdit, 'String'));
guidata(hObject, data);

% --- Executes during object creation, after setting all properties.
function CDCEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%function ETTEdit_Callback(hObject, eventdata, handles)
%data = guidata(hObject);
%data.ETTable = str2num(get(handles.ETTEdit, 'String'));
%guidata(hObject, data);

% --- Executes during object creation, after setting all properties.
%function ETTEdit_CreateFcn(hObject, eventdata, handles)
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%    set(hObject,'BackgroundColor','white');
%end


% --- Executes on button press in PlantLoadButton.
function PlantLoadButton_Callback(hObject, eventdata, handles)
data = guidata(hObject);
filename = inputdlg("Specify the file name for plant model", "Open file");
load(filename{1,1});
data.task_freq = freq;
data.A = A;
data.B2= B;             %control input matrix
data.B1= Bd;            %disturbance matrix
data.Q = Q;
data.R = R;
%num_disturbance = size(Bd,2);
num_output      = size(B,2);
data.n_tasks = num_output;
data.reserved_util = 0.;
%data.num_output = num_output;
%data.num_disturbance = num_disturbance;
guidata(hObject, data);


% --- Executes on button press in LowButton.
function LowButton_Callback(hObject, eventdata, handles)
data = guidata(hObject);
data.Threshold = 50;
data.PerText1 = 'Low';
guidata(hObject, data);
set(handles.DevVal, 'String', num2str(data.Threshold));


% --- Executes on button press in MediumButton.
function MediumButton_Callback(hObject, eventdata, handles)
data = guidata(hObject);
data.Threshold = 25;
data.PerText1 = 'Medium';
guidata(hObject, data);
set(handles.DevVal, 'String', num2str(data.Threshold));


% --- Executes on button press in HighButton.
function HighButton_Callback(hObject, eventdata, handles)
data = guidata(hObject);
data.Threshold = 5;
data.PerText1 = 'High';
guidata(hObject, data);
set(handles.DevVal, 'String', num2str(data.Threshold));


% --- Executes on button press in TuneButton.
function TuneButton_Callback(hObject, eventdata, handles)
data = guidata(hObject);
disp('<< == Start tuning the K matrix! == >>');
A = data.A;
B1 = data.B1;
B2 = data.B2;
%Q = eye(size(A));
%num_output = data.num_output;
%R = eye(num_output,num_output);
Q = data.Q;
R = data.R;
task_freq = data.task_freq;
proc_lim = data.proc_lim;
reserved_util = data.reserved_util;
n_tasks = data.n_tasks;

max_iter = 700;

fprintf('%s\t\t%s\t\t%s\t%s\n', 'Iter', 'nnz', 'H2_norm', 'Est Util');
gamma = 0;
%% A larger optimization rate to reduce iterations
%gamma_inc = 0.002;
%% A smaller optimization rate for fine-tuning
gamma_inc = 1e-04;

fprintf('%d\t\t', 0);
% dense optimal prototype
soln_dense_opt = H2sparse(A,B1,B2,Q,R,gamma);
soln = soln_dense_opt;
K = soln.F;

while (1)
iter = 1;
while (true)
    % Demultiplexor, adopt as necessary
    blockified_tasks = Demultiplexor(K, n_tasks);

    % ET_Estimator, takes in reshaped K matrix (1*x) (blockified_tasks)
    ET = ET_Estimator(blockified_tasks, data.PerfMat);

    % RMA Consolidator, Feed Back: tasks should be simplified a total factor
    % of shrink_factor (in terms of ET), to be fit into the computation nodes.
    [procs_occupied, tasks_on_proc, proc_util, proc_alloc_to_task, unfulfilled_tasks, ~] = ...
        RMA_Consolidator(proc_lim, 1. - reserved_util, ET, task_freq);
    
    util = sum(ET .* task_freq)/procs_occupied;
    fprintf ('%f\n',util);
    
    gamma_g(iter) = gamma;
    nnz_g(iter) = soln.nnz;
    deviation_g(iter) = abs(soln.J - soln_dense_opt.J) / soln_dense_opt.J;
    util_g(iter) = util;

    if (~isempty(unfulfilled_tasks))
        % needs to further simplify
        gamma = gamma+gamma_inc;
        fprintf('%d\t\t', iter);
        soln = H2sparse(A,B1,B2,Q,R,gamma);
        K = soln.F;
        % Check iteration boundary
        if(iter > max_iter)
			error 'Error: Iteration exhausted. Hardware constraint not met!';
        end
        iter = iter + 1;
    else
        break;
    end
end

disp('Verifying Solution...');
RMA_Verify(procs_occupied, tasks_on_proc, proc_util, proc_alloc_to_task, ET, 1. - reserved_util, task_freq);

deviation = abs(soln.J - soln_dense_opt.J) / soln_dense_opt.J * 100;
set(handles.NOCUEdit, 'String', procs_occupied);
set(handles.DEVEdit, 'String', deviation);
set(handles.AvgUtil, 'String', util*100);

    gamma_g(iter) = gamma;
    nnz_g(iter) = soln.nnz;
    deviation_g(iter) = deviation/100;
    util_g(iter) = util;
    
    figure;
    plot(gamma_g,deviation_g*100,'DisplayName','deviation_g');
    hold on;
    plot(gamma_g,nnz_g,'DisplayName','nnz_g');
    hold off;

if deviation > data.Threshold
    % Deviation requirement not met, pop up selection prompt.
    set(handles.DEVEdit, 'ForegroundColor', [1.0, 0.0, 0.0]);
    % Construct a questdlg with two options
    choice = questdlg('WARNING: Deviation constraint not met! You can ...', ...
	'Alternatives', ...
	'Increase a Processing Node','Continue with Higher Deviation','Increase a Processing Node');
    % Handle response
    switch choice
        case 'Increase a Processing Node'
            data.proc_lim = data.proc_lim+1;
            proc_lim = data.proc_lim;
            set(handles.NOACEdit, 'String', data.proc_lim);
            loop = 0;
            gamma = gamma_g(find(util_g <= proc_lim/(proc_lim-1), 1) - 1);
            gamma_inc = gamma_inc * 0.5;
            disp('<< == Restart tuning the K matrix == >>');
            fprintf('%d\t\t', loop);
            soln = H2sparse(A,B1,B2,Q,R,gamma);
            K = soln.F;
        case 'Continue with Higher Deviation'
            break;
    end
else
    set(handles.DEVEdit, 'ForegroundColor', [0.0, 0.0, 0.0]);
    msgbox('Both hardware and deviation constraints are met.');
    break;
end
end

% Packing Data for GUI.
data.PercentDeviation = deviation;
data.procs_occupied = procs_occupied;
data.tasks_on_proc = tasks_on_proc;
data.AvgUtil = util;
data.proc_alloc_to_task = proc_alloc_to_task;
data.blockified_tasks = blockified_tasks;

data.utilNode = proc_util;
data.Kref = soln_dense_opt.F;
data.K1 = soln.F;

guidata(hObject, data);

% --- Executes on button press in StopButton.
function StopButton_Callback(hObject, eventdata, handles)
disp('<< == Finish the sparsity routine! == >>');
clear all;
close all;
clc;
close;



function NOCUEdit_Callback(hObject, eventdata, handles)
data = guidata(hObject);
set(handles.NOCUEdit, 'String', data.procs_occupied);

% --- Executes during object creation, after setting all properties.
function NOCUEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function DEVEdit_Callback(hObject, eventdata, handles)
data = guidata(hObject);
set(handles.DEVEdit, 'String', num2str(data.PercentDeviation));

% --- Executes during object creation, after setting all properties.
function DEVEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RunOnHW.
function RunOnHW_Callback(hObject, eventdata, handles)
% hObject    handle to RunOnHW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = guidata(hObject);
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
procs_occupied = data.procs_occupied;
n_tasks = data.n_tasks;
tasks_on_proc = data.tasks_on_proc;
proc_alloc_to_task = data.proc_alloc_to_task;
blockified_tasks = data.blockified_tasks;

PubPort = 60001:(60000+max_num_of_processors);
TaskPort = 50001:(50000+n_tasks);
bufSize = uint32([0 0]);
for i = 1:procs_occupied
    PubComm(i) = tcpip('0.0.0.0',PubPort(i),...
                       'NetworkRole','server',...
                       'ByteOrder','littleEndian',...
                       'Timeout',10.0);
end
disp ('Waiting for Connection...');
for i = 1:procs_occupied
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
for i = 1:procs_occupied
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
guidata(hObject, data);


function DevVal_Callback(hObject, eventdata, handles)
data = guidata(hObject);
data.Threshold = str2num(get(handles.DevVal, 'String'));
guidata(hObject, data);


% --- Executes during object creation, after setting all properties.
function DevVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function AvgUtil_Callback(hObject, eventdata, handles)
data = guidata(hObject);
set(handles.AvgUtil, 'String', data.AvgUtil);


% --- Executes during object creation, after setting all properties.
function AvgUtil_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ViewResp.
function ViewResp_Callback(hObject, eventdata, handles)
% hObject    handle to ViewResp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = guidata(hObject);
noise = rand(size(data.B1,2), 1000);
respref = zeros(size(data.A,1), 1000);
resp = zeros(size(data.A,1), 1001);
for t = 1:1000
    xdref = data.A*respref(:,t) + data.B1*noise(:,t) - data.B2*data.Kref(1:4,:)*(respref(:,t)-[zeros(41,1);1;zeros(39,1);1]);
    respref(:,t+1) = respref(:,t)+xdref*0.025;
    xd = data.A*resp(:,t) + data.B1*noise(:,t) - data.B2*data.K1(1:4,:)*(resp(:,t)-[zeros(41,1);1;zeros(39,1);1]);
    resp(:,t+1) = resp(:,t)+xd*0.025;
end
figure;
plot([0:0.025:25], resp([40,42,81,82],:)','LineWidth',1.5); hold on; plot([0:0.025:25], respref([40,42,81,82],:)'); hold off;
grid on;

% --- Executes on button press in ViewHWUtil.
function ViewHWUtil_Callback(hObject, eventdata, handles)
% hObject    handle to ViewHWUtil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = guidata(hObject);
figure; bar(data.utilNode);

% --- Executes on button press in LoadHWPerf.
function LoadHWPerf_Callback(hObject, eventdata, handles)
% hObject    handle to LoadHWPerf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    data = guidata(hObject);
    filename = inputdlg("Specify the file name for hardware capacity profile", "Open file");
    load(filename{1,1}, 'wcet_table');
    data.PerfMat = wcet_table;
    guidata(hObject, data);
