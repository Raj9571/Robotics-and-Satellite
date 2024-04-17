
X = 6; %
p1 = 3.31 + X/30;
p2 = 0.116 + X/500;
p3 = 0.16 + X/400;
m1 = 2.5;
m2 = 1.5;
l = 0.4;
g = 9.81;
qr = [0.4; 0];
% Declare global variable to store control input
global control_input_history;
control_input_history = [];

% Define suitable D and Kp matrices and initial conditions
D = diag([0.1, 0.1]); % Example values for damping matrix
Kp = diag([1, 1]); 
q0 = [1; -1]; % Initial joint angles
qd0 = [1; 1]; % Initial joint velocities


e0 = q0 - qr;
edot0 = qd0; 

initial_conditions = [e0; edot0];
time_span = [0 50];


dynamics = @(t, y) manipulator_error_dynamics(t, y, p1, p2, p3, m1, m2, l, g, D, Kp, qr);


[t, state] = ode45(dynamics, time_span, initial_conditions);


control_input_history = zeros(length(t), 2); 

% Compute control inputs using the state history
for i = 1:length(t)
    [~, control_input] = manipulator_error_dynamics(t(i), state(i, :)', p1, p2, p3, m1, m2, l, g, D, Kp, qr);
    control_input_history(i, :) = control_input';
end

% Extract error and error derivative from the state
e = state(:, 1:2);    % Error in position
edot = state(:, 3:4); % Error in velocity

% actual q and qdot from error and error derivative
q = e + qr';
qd = edot;

%  joint angles
figure;
plot(t, q);
title('Joint Angles');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('q1', 'q2');

%  joint velocities
figure;
plot(t, qd);
title('Joint Velocities');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('qd1', 'qd2');

%  control inputs
figure;
plot(t, control_input_history);
title('Control Inputs');
xlabel('Time (s)');
ylabel('Control Input (u)');
legend('u1', 'u2');

% Dynamics function
function [dedt, u] = manipulator_error_dynamics(t, y, p1, p2, p3, m1, m2, l, g, D, Kp, qr)
    e1 = y(1:2);
    e2 = y(3:4);
    q = e1 + qr;
    q_dot = e2;
    M = inertia_matrix(q, p1, p2, p3);
    C = coriolis_matrix(q, q_dot, p3);
    g_q = gravity_effect(q, m1, m2, l, g);
    nu = - tanh(e2);
    u = g_q - Kp * e1 + nu;
    e_dot_dot = M \ (u - C * e2 - D * e2 - g_q);
    dedt = [e2; e_dot_dot];
end


% Inertia matrix function
function M = inertia_matrix(q, p1, p2, p3)
    M = [p1 + 2*p3*cos(q(2)), p2 + p3*cos(q(2));
         p2 + p3*cos(q(2)), p2];
end

% Coriolis matrix function
function C = coriolis_matrix(q, qd, p3)
    C = p3*sin(q(2)) * [-qd(2), qd(1) + qd(2);
                         qd(1), 0];
end

% Gravity effect function
function g_q = gravity_effect(q, m1, m2, l, g)
    g_q = g * [l*(m1 + m2)*cos(q(1)) + l*m2*cos(q(1) + q(2));
               l*m2*cos(q(2))];
end
