% test_flight_dynamics.m
clear all; close all; clc;

% Load aircraft configuration
aircraft = init_aircraft();

% Time settings
dt = 0.1;
t_final = 20; % 20 seconds total
time = 0:dt:t_final;
N = length(time);

% Initialize storage
x = zeros(1,N); y = zeros(1,N); z = zeros(1,N);
V = zeros(1,N); L = zeros(1,N); D = zeros(1,N);
weight = aircraft.mass * 9.81; % For comparison

% Initial conditions - matching flight_dynamics.m requirements
current_state.x = 0;
current_state.y = 0;
current_state.z = 1000; % Start at 1000m
current_state.V = 15;   % Initial velocity
current_state.gamma = 0; % Level flight initially
current_state.phi = 0;   % No initial bank
current_state.psi = 0;   % Initial heading

% Store initial conditions
x(1) = current_state.x;
y(1) = current_state.y;
z(1) = current_state.z;
V(1) = current_state.V;

% Wind (zero for basic testing)
wind.Wx = 0; wind.Wy = 0; wind.U = 0;

% Test two phases:
% 1. Powered level flight (0-10s)
% 2. Unpowered glide (10-20s)
for i = 2:N
    t = time(i);
    
    % Control inputs matching flight_dynamics.m requirements
    if t < 10
        control_input.thrust = 5; % Some thrust for level flight
    else
        control_input.thrust = 0; % No thrust for gliding
    end
    control_input.phi_command = 0;    % Maintain level bank
    
    % Simulate and store results
    next_state = flight_dynamics(current_state, control_input, dt, wind, aircraft);
    x(i) = next_state.x;
    y(i) = next_state.y;
    z(i) = next_state.z;
    V(i) = next_state.V;
    L(i) = next_state.debug.L;
    D(i) = next_state.debug.D;
    current_state = next_state;
end

% Plot results
figure(1)
subplot(2,2,1)
plot(time, z)
ylabel('Altitude (m)')
title('Altitude vs Time')
grid on
hold on
plot([10 10], ylim, 'r--') % Mark thrust cutoff
legend('Flight Path', 'Thrust Cutoff')

subplot(2,2,2)
plot(time, V)
ylabel('Airspeed (m/s)')
title('Airspeed vs Time')
grid on
hold on
plot([10 10], ylim, 'r--')

subplot(2,2,3)
plot(x, z)
xlabel('Distance (m)')
ylabel('Altitude (m)')
title('Flight Path')
grid on
axis equal

subplot(2,2,4)
plot(time, L, 'b-', time, D, 'r-', time, ones(size(time))*weight, 'k--')
ylabel('Forces (N)')
title('Forces vs Time')
legend('Lift', 'Drag', 'Weight')
grid on
hold on
plot([10 10], ylim, 'r--')

% Print analysis
fprintf('\nFlight Analysis:\n')

% Powered phase
t_idx = time <= 10;
fprintf('\nPowered Flight (0-10s):\n')
fprintf('Average Speed: %.1f m/s\n', mean(V(t_idx)))
fprintf('Altitude Change: %.1f m\n', z(find(t_idx,1,'last')) - z(1))
fprintf('Average L/W ratio: %.2f\n', mean(L(t_idx))/weight)
fprintf('Average L/D ratio: %.1f\n', mean(L(t_idx))/mean(D(t_idx)))

% Gliding phase
t_idx = time > 10;
fprintf('\nGliding Flight (10-20s):\n')
fprintf('Average Speed: %.1f m/s\n', mean(V(t_idx)))
fprintf('Altitude Change: %.1f m\n', z(end) - z(find(t_idx,1)))
fprintf('Average L/W ratio: %.2f\n', mean(L(t_idx))/weight)
fprintf('Average L/D ratio: %.1f\n', mean(L(t_idx))/mean(D(t_idx)))

% Calculate glide ratio
dx = x(end) - x(find(t_idx,1));
dz = z(end) - z(find(t_idx,1));
fprintf('Glide ratio: %.1f:1\n', abs(dx/dz))