% test_flight_dynamics_2.m
clear all; close all; clc;

% Load aircraft configuration
aircraft = init_aircraft();

% Time settings
dt = 0.1;
t_final = 100; % 30 seconds for each test
time = 0:dt:t_final;
N = length(time);

% Initialize arrays for 3D tracking
x = zeros(1,N); y = zeros(1,N); z = zeros(1,N);
V = zeros(1,N); psi = zeros(1,N); phi = zeros(1,N);
gamma = zeros(1,N);

% Wind (zero for testing)
wind.Wx = 0; wind.Wy = 0; wind.U = 0;

fprintf('\nFlight Dynamics Test Based on Glider Handbook:\n')

%% Test 1: Best Glide
current_state.x = 0;
current_state.y = 0;
current_state.z = 1000;            % Starting altitude
current_state.V = aircraft.cruise_speed;  % Use aircraft's cruise speed (15 m/s)
current_state.gamma = -0.025;       % Slight negative angle for glide (-3 degrees)
current_state.phi = 0;             % No bank
current_state.psi = 0;             % Initial heading

control_input.thrust = 0;          % No thrust for gliding
control_input.phi_command = 0;     % No banking command
control_input.gamma_command = -0.025; % Match initial gamma for steady glide

[x1,y1,z1,V1,psi1] = run_test_3d(current_state, control_input, time, dt, wind, aircraft);
analyse_glide(x1, z1, V1, time);

% Plot gliding test
figure(1)
plot3D_trajectory(x1, y1, z1, V1, aircraft, 'Gliding Test');

%% Test 2: Turning Flight
fprintf('\nTest 2: Turning Flight Performance\n')
current_state.x = 0;
current_state.y = 0;
current_state.z = 1000;
current_state.phi = 0;
current_state.psi = -0.05;
control_input.phi_command = 30 * pi/180; % 30 degree bank based on literature
control_input.gamma_command = -0.05;

[x2,y2,z2,V2,phi2] = run_turn_test_3d(current_state, control_input, time, dt, wind, aircraft);
analyse_turn(phi2, V2, time);

% Plot turning test
figure(2)
plot3D_trajectory(x2, y2, z2, V2, aircraft, 'Turning Test');

%% Helper Functions
function [x,y,z,V,psi] = run_test_3d(init_state, control, time, dt, wind, aircraft)
    N = length(time);
    x = zeros(1,N); y = zeros(1,N); z = zeros(1,N); 
    V = zeros(1,N); psi = zeros(1,N);
    current_state = init_state;
    
    for i = 1:N
        x(i) = current_state.x;
        y(i) = current_state.y;
        z(i) = current_state.z;
        V(i) = current_state.V;
        psi(i) = current_state.psi;
        current_state = flight_dynamics(current_state, control, dt, wind, aircraft);
    end
end

function [x,y,z,V,phi] = run_turn_test_3d(init_state, control, time, dt, wind, aircraft)
    N = length(time);
    x = zeros(1,N); y = zeros(1,N); z = zeros(1,N); 
    V = zeros(1,N); phi = zeros(1,N);
    current_state = init_state;
    
    for i = 1:N
        x(i) = current_state.x;
        y(i) = current_state.y;
        z(i) = current_state.z;
        V(i) = current_state.V;
        phi(i) = current_state.phi;
        current_state = flight_dynamics(current_state, control, dt, wind, aircraft);
    end
end

function plot3D_trajectory(x, y, z, V, aircraft, titleStr)

    %% Set up global plotting settings
    set(0, 'DefaultAxesFontSize', 12);
    set(0, 'DefaultFigureColor', 'w');
    set(0, 'DefaultTextInterpreter', 'tex');
    set(0, 'DefaultAxesFontName', 'Times New Roman');
    set(0, 'DefaultTextFontName', 'Times New Roman');

    % Create figure with thesis-appropriate dimensions
    fig_handle = figure;
    set(fig_handle, 'Units', 'centimeters');
    set(fig_handle, 'Position', [0 0 17 10]); % Width and height in centimeters
    hold on
    
    % Create velocity-based coloring
    normalized_velocities = (V - aircraft.min_speed) / ...
        (aircraft.max_speed - aircraft.min_speed);
    colors = jet(100); % Create colormap
    
    % Plot path segments with velocity colors
    for i = 1:length(V)-1
        color_idx = max(1, min(100, round(normalized_velocities(i) * 99) + 1));
        plot3([x(i), x(i+1)], [y(i), y(i+1)], [z(i), z(i+1)], ...
            'Color', colors(color_idx,:), 'LineWidth', 2);
    end
    
    % Plot start and end points (plot these after the path)
    h_start = plot3(x(1), y(1), z(1), 'go', 'MarkerSize', 10, 'LineWidth', 2); % Start point
    h_end = plot3(x(end), y(end), z(end), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % End point
    
    % Add colorbar
    c = colorbar;
    c.Label.String = 'Velocity (m/s)';
    clim([aircraft.min_speed aircraft.max_speed]);
    
    % Set axis limits with some padding
    x_range = max(x) - min(x);
    y_range = max(y) - min(y);
    z_range = max(z) - min(z);
    
    x_padding = max(x_range * 0.1, 1);
    y_padding = max(y_range * 0.1, 1);
    z_padding = max(z_range * 0.1, 1);
    
    xlim([min(x) - x_padding, max(x) + x_padding]);
    ylim([min(y) - y_padding, max(y) + y_padding]);
    zlim([min(z) - z_padding, max(z) + z_padding]);
    
    % Add labels and title
    xlabel('X Distance (m)');
    ylabel('Y Distance (m)');
    zlabel('Altitude (m)');
    title(titleStr);
    grid on
    axis equal
    view(45, 30)
    
    % Add legend with the point markers
    legend([h_start, h_end], {'Start', 'End'}, 'Location', 'best');
    
    hold off
end


function analyse_glide(x, z, V, time)
    % Original performance metrics
    dx = x(end) - x(1);
    dz = z(1) - z(end);
    glide_ratio = abs(dx/dz);
    avg_speed = mean(V);
    sink_rate = dz/time(end);
    
    fprintf('Glide Analysis:\n');
    fprintf('Initial Speed: %.1f m/s\n', V(1));
    fprintf('Final Speed: %.1f m/s\n', V(end));
    fprintf('Average Speed: %.1f m/s\n', avg_speed);
    fprintf('Min Speed: %.1f m/s\n', min(V));
    fprintf('Max Speed: %.1f m/s\n', max(V));
    fprintf('Sink Rate: %.2f m/s\n', sink_rate);
    fprintf('Glide Ratio: %.1f:1\n', glide_ratio);
    fprintf('Total Altitude Lost: %.1f m\n', dz);
    
    % Print velocity progression every second
    fprintf('\nVelocity Progression (every second):\n');
    fprintf('Time (s)\tVelocity (m/s)\n');
    samples_per_second = round(1/time(2)); % Assuming uniform time steps
    for i = 1:samples_per_second:length(time)
        fprintf('%.1f\t\t%.1f\n', time(i), V(i));
    end
end

function analyse_turn(phi, V, time)
    avg_bank = mean(phi(end-100:end));  % Average over last 10 seconds
    load_factor = 1/cos(avg_bank);
    turn_radius = mean(V(end-100:end))^2 / (9.81 * tan(avg_bank));
    
    fprintf('Turn Analysis:\n');
    fprintf('Bank Angle: %.1f degrees\n', avg_bank * 180/pi);
    fprintf('Load Factor: %.2f G\n', load_factor);
    fprintf('Turn Radius: %.1f m\n', turn_radius);
end

