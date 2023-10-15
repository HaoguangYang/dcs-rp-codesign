% Platoon controller: state space = [
%     lengths to center plate;
%     heading error;
%     base_local_x,
%     base_local_y,
%     plate heading angle
% ]

dt = 1/50;
% initial state to linearize the problem
l10 = 1.0;
l20 = 1.0;
l30 = 1.0;

delta = pi*2/3;
A = eye(9);
B = [sqrt(l10^2+(dt/3)^2)-l10, -dt, 0, sin(delta)*dt/2, 0, 0, -sin(delta)*dt/2, 0, 0;
    -sin(delta)*dt/2, 0, 0, sqrt(l20^2+(dt/3)^2)-l20, -dt, 0, sin(delta)*dt/2, 0, 0;
    sin(delta)*dt/2, 0, 0, -sin(delta)*dt/2, 0, 0, sqrt(l30^2+(dt/3)^2)-l30, -dt, 0;
    -atan(dt/(3*l10)), 0, dt, -atan(dt/(3*l20)), 0, 0, -atan(dt/(3*l30)), 0, 0;
    -atan(dt/(3*l10)), 0, 0, -atan(dt/(3*l20)), 0, dt, -atan(dt/(3*l30)), 0, 0;
    -atan(dt/(3*l10)), 0, 0, -atan(dt/(3*l20)), 0, 0, -atan(dt/(3*l30)), 0, dt;
    1.0*dt*2/3, 0, 0, cos(delta)*dt*2/3, 0, 0, cos(-delta)*dt*2/3, 0, 0;
    0, 0, 0, sin(delta)*dt*2/3, 0, 0, sin(-delta)*dt*2/3, 0, 0;
    atan(dt/(3*l10)), 0, 0, atan(dt/(3*l20)), 0, 0, atan(dt/(3*l30)), 0, 0];

Q = diag([500, 500, 500, 100, 100, 100, 500, 500, 100]);
R = 10*eye(9);
K = dlqr(A, B, Q, R);

cont_sys = d2c(ss(A, B, eye(9), zeros(9,9), dt));
% This continuous system matrix can be used to design controller.
A_cont = cont_sys.A;
B_cont = cont_sys.B;

A = A_cont;
B = B_cont;
Bd = 0.2*B;

tau = dt * ones(9,1);

save("sys_coordinated_robots.mat", "A", "B", "Bd", "Q", "R", "tau");
