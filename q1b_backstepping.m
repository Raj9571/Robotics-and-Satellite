% Define the parameters
X = 6;
p1 = 3.31 + X/30;
p2 = 0.116 + X/500;
p3 = 0.16 + X/400;

qr = 0.3; % Desired reference position

% Initial conditions
q0 = 0; 
q_dot0 = 0; 
xi0 = 0;
% Control law parameters
k1 = 15; % Control gain for e1
e0=q0-qr;
% Time span for the simulation
tspan = [0, 50]; % 10 seconds simulation


initial_conditions = [e0; q_dot0; xi0];

% Anonymous function to pass the extra parameters to the ODE function
system_dynamics = @(t, y) system_dynamics_impl(t, y, p1, p2, p3, qr, k1);

% Simulate the system using the anonymous function
[t, states] = ode45(system_dynamics, tspan, initial_conditions);

% Extract the joint angle, velocity, and xi from the states
q = states(:, 1)+qr;
q_dot = states(:, 2);
xi = states(:, 3);


tau_history = zeros(length(t), 1);
for i = 1:length(t)
    [~, ~, tau_history(i)] = system_dynamics_impl(t(i), states(i, :), p1, p2, p3, qr, k1);
end

% Plot the joint angle over time
figure;
plot(t, q, 'b', t, q_dot, 'r', t, qr*ones(size(t)), 'r--');
xlabel('Time (s)');
ylabel('Joint angle (rad) and velocity (rad/s)');
title('Joint Angle and Velocity Tracking');
legend('Joint angle q', 'Joint velocity q\_dot', 'Desired joint angle q\_r');

% Plot the control input over time
figure;
plot(t, tau_history, 'k-');
xlabel('Time (s)');
ylabel('Control input tau (N*m)');
title('Control Input History');

% Define the system dynamics function
function [dydt, xi, tau] = system_dynamics_impl(t, e, p1, p2, p3, qr, k1)
    % Extract the current state
     q = e(1)+qr;
    q_dot = e(2);
    xi = e(3);

  
    m_q = p1 + 2*p3*cos(q);
    c_q_q_dot = -p3*sin(q)*q_dot;

    %  control law here based on the backstepping design
  tau =-e(2)-e(2)^2*sin(q)*p3-p3*sin(q)*(xi+e(1)*e(2)*p3*sin(q))/m_q-p3*e(1)*e(2)^2*cos(q)-(xi+e(1)*e(2)*p3*sin(q))/m_q-k1*(xi+e(1)+e(2)+p3*e(1)*e(2)*sin(q)) ; % Actual control input
  
    q_ddot = (xi - c_q_q_dot * q_dot) / m_q; 
    xi_dot = tau; 
    % Return the state derivatives and control input
    dydt = [q_dot; q_ddot; xi_dot];
  
end
