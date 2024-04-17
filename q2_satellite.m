
    % Parameters
    alpha1 =  (126)/100; 
    J = [20 + alpha1, 1.2, 0.9; 1.2, 17 + alpha1, 1.4; 0.9, 1.4, 15 + alpha1];
    k1 = 2; 
    k2 = 2;
    % Initial conditions
    rho_0 = [-0.02, 0, 0.045]';
    omega_0 = [0.004, -0.007, 0.017]';
    y0 = [rho_0; omega_0];

   
    tspan = [0, 50]; 

    [T, Y] = ode45(@(t, y) rigid_body_ode(t, y, J, k1, k2), tspan, y0);

    % Extract rho and omega from Y
    rho = Y(:, 1:3);
    omega = Y(:, 4:6);

    % Control input u
    nu=-k1*tanh(omega);
    u=nu - 2 * k2 * rho;
   
    % Plot the results
    figure;
    subplot(2,1,1);
    plot(T, rho);
    title('Time history of rho');
    xlabel('Time (s)');
    ylabel('rho');
    legend('rho_1', 'rho_2', 'rho_3');

    subplot(2,1,2);
    plot(T, omega);
    title('Time history of omega');
    xlabel('Time (s)');
    ylabel('omega');
    legend('omega_1', 'omega_2', 'omega_3');

    figure;
    plot(T, u);
    title('Control input u');
    xlabel('Time (s)');
    ylabel('u');
    legend('u_1', 'u_2', 'u_3');


function dydt = rigid_body_ode(t, y, J, k1, k2)
    
    rho = y(1:3);
    omega = y(4:6);

    
    I = eye(3);
    nu=-k1*tanh(omega);
    u=nu - 2 * k2 * rho;
    % Equations of motion
    rho_dot = (I + skew_symmetric(rho) + rho*rho') * omega;
    omega_dot = J \  (-skew_symmetric(omega) *J*omega +u );

    % Concatenate derivatives
    dydt = [rho_dot; omega_dot];
end

function S = skew_symmetric(v)
    % Returns the skew-symmetric matrix of a 3x1 vector
    S = [  0   -v(3)  v(2);
          v(3)  0   -v(1);
         -v(2)  v(1)  0];
end
