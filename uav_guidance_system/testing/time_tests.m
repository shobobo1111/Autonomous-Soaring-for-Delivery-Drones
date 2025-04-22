% % test_thermal_timing_comparison.m
% % This script systematically tests the effect of varying delivery time 
% % constraints on path planning in thermal conditions
% 
% clear all; close all; clc;
% 
% % Check if init_mission.m accepts parameters
% try
%     test_mission = init_mission(250);
%     if test_mission.time_to_deliver ~= 250
%         error('init_mission.m does not properly handle input parameters. Please update it according to the instructions.');
%     end
%     fprintf('init_mission.m successfully accepts time parameters.\n\n');
% catch e
%     error(['Error: init_mission.m needs to be modified to accept a time parameter.\n', ...
%            'Please update it to have the following structure:\n', ...
%            'function mission = init_mission(varargin)\n', ...
%            '    mission.time_to_deliver = 300;\n', ...
%            '    mission.alt_floor = 150;\n', ...
%            '    mission.alt_floor_tol = 50;\n', ...
%            '    if nargin > 0 && ~isempty(varargin{1})\n', ...
%            '        mission.time_to_deliver = varargin{1};\n', ...
%            '    end\n', ...
%            'end\n']);
% end
% 
% % Load configurations
% aircraft = init_aircraft();
% battery = init_battery();
% waypoint = init_waypoint();
% 
% %% Create wind field
% % Define grid to match flight space
% x_limits = [-100, 2000];
% y_limits = [-500, 500];
% z_limits = [0, 800];
% grid_spacing = 20;
% 
% % Create wind field structure
% wind_field = init_wind_field(x_limits, y_limits, z_limits, grid_spacing);
% 
% % Add thermal pattern - using 'simple_sys' for consistent thermal placement
% wind_field = create_thermals(wind_field, 'simple_sys');
% 
% %% Initial conditions
% initial_state.x = 0;
% initial_state.y = 0;
% initial_state.z = 300;
% initial_state.V = 15;        % m/s
% initial_state.V_rel = 15;    % m/s
% initial_state.gamma = 0;     % rad
% initial_state.phi = 0;       % rad
% initial_state.psi = 0;       % rad
% 
% %% Set up test parameters
% % Define range of delivery times to test
% min_delivery_time = 100;  % seconds
% max_delivery_time = 400;  % seconds
% num_tests = 2;
% 
% delivery_times = linspace(min_delivery_time, max_delivery_time, num_tests);
% num_steps = 10;  % Maximum steps per simulation
% 
% % Storage for all test results
% all_results = cell(num_tests, 1);
% 
% % Create a special colormap for time visualization (blue to red gradient)
% % Using viridis-like colormap for better visualization 
% [custom_cmap, ~] = init_colour_viridis();
% 
% %% Run all tests
% for test_idx = 1:num_tests
%     fprintf('\n\n==================================================\n');
%     fprintf('STARTING TEST %d of %d: Time to deliver = %.1f seconds\n', ...
%         test_idx, num_tests, delivery_times(test_idx));
%     fprintf('==================================================\n\n');
% 
%     % Run simulation with current time parameter
%     complete_path = generate_multi_step_path(initial_state, aircraft, wind_field, num_steps, battery, delivery_times(test_idx));
% 
%     % Store results
%     all_results{test_idx} = complete_path;
%     all_results{test_idx}.time_to_deliver = delivery_times(test_idx);
% 
%     % Calculate summary statistics for this run
%     total_flight_time = complete_path.steps * 2; % seconds
%     max_altitude = max(complete_path.states.z);
%     final_distance = sqrt((waypoint.x - complete_path.states.x(end))^2 + ...
%                          (waypoint.y - complete_path.states.y(end))^2);
% 
%     % Calculate final battery percentage
%     battery_capacity = battery.capacity_Wh;
%     base_power = battery.base_power;
%     thrust_coefficient = battery.thrust_coefficient;
% 
%     % Simple battery calculation to estimate final percentage
%     total_energy_used = 0;
%     for s = 1:length(complete_path.states.thrust)
%         thrust = complete_path.states.thrust(s);
%         power = base_power + thrust_coefficient * thrust;
%         energy = power * (2/3600); % 2 seconds converted to hours
%         total_energy_used = total_energy_used + energy;
%     end
% 
%     remaining_battery = battery_capacity - total_energy_used;
%     battery_percentage = (remaining_battery / battery_capacity) * 100;
% 
%     all_results{test_idx}.final_battery = remaining_battery;
%     all_results{test_idx}.final_battery_percentage = battery_percentage;
% 
%     fprintf('\nTest %d Summary:\n', test_idx);
%     fprintf('  Time constraint: %.1f seconds\n', delivery_times(test_idx));
%     fprintf('  Actual flight time: %.1f seconds\n', total_flight_time);
%     fprintf('  Max altitude reached: %.1f meters\n', max_altitude);
%     fprintf('  Final distance to waypoint: %.1f meters\n', final_distance);
%     fprintf('  Remaining battery: %.2f%% (%.2f Wh out of %.2f Wh)\n', battery_percentage, remaining_battery, battery_capacity);
%     fprintf('  Steps completed: %d\n', complete_path.steps);
% end
% 
% %% Set up global plotting settings
% set(0, 'DefaultAxesFontSize', 12);
% set(0, 'DefaultFigureColor', 'w');
% set(0, 'DefaultTextInterpreter', 'tex');
% set(0, 'DefaultAxesFontName', 'Times New Roman');
% set(0, 'DefaultTextFontName', 'Times New Roman');
% 
% %% Create visualization of all paths
% % Create a figure for the 3D paths with thesis-appropriate dimensions
% fig_handle = figure;
% set(fig_handle, 'Units', 'centimeters');
% set(fig_handle, 'Position', [0 0 17 13]); % Width and height in centimeters
% 
% % Plot the wind field first (light gray)
% hold on;
% downsample = 3;  % Adjust for density of wind vectors
% vertical_threshold = 0.2;  % Only show significant vertical winds
% significant_vertical = abs(wind_field.U) >= vertical_threshold;
% 
% % Create a combined mask that both downsamples and checks for significance
% mask = false(size(wind_field.X));
% mask(1:downsample:end, 1:downsample:end, 1:downsample:end) = true;
% mask = mask & significant_vertical;
% 
% % Plot wind vectors
% quiver3(wind_field.X(mask), wind_field.Y(mask), wind_field.Z(mask), ...
%         wind_field.Wx(mask), wind_field.Wy(mask), wind_field.U(mask), ...
%         'Color', [0.8 0.8 0.8], 'LineWidth', 0.8, 'MaxHeadSize', 0.3);
% 
% % Plot waypoint marker
% plot3(waypoint.x, waypoint.y, waypoint.z, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
% 
% % Plot all paths with color coding for time to deliver
% cmap_indices = round(linspace(1, size(custom_cmap, 1), num_tests));
% 
% for test_idx = 1:num_tests
%     path = all_results{test_idx};
% 
%     % Get the color for this path based on time constraint
%     path_color = custom_cmap(cmap_indices(test_idx), :);
% 
%     % Plot this path
%     plot3(path.states.x, path.states.y, path.states.z, ...
%         'Color', path_color, 'LineWidth', 1.5);
% 
%     % Mark the start point for the first path only
%     if test_idx == 1
%         plot3(initial_state.x, initial_state.y, initial_state.z, ...
%             'go', 'MarkerSize', 10, 'LineWidth', 2);
%     end
% 
%     % Mark the end point of each path with a dot
%     plot3(path.states.x(end), path.states.y(end), path.states.z(end), ...
%         'o', 'Color', path_color, 'MarkerSize', 8, 'LineWidth', 1.5);
% end
% 
% % Format plot
% xlabel('East (m)', 'FontWeight', 'bold');
% ylabel('North (m)', 'FontWeight', 'bold');
% zlabel('Altitude (m)', 'FontWeight', 'bold');
% title('Effect of Delivery Time Constraints on Path Planning', 'FontWeight', 'bold');
% 
% % Set position of axes in the figure for consistency
% set(gca, 'Units', 'centimeters');
% set(gca, 'Position', [2 1.5 11 10]);  % [left bottom width height]
% 
% % Set view and grid
% grid off;
% view(30, 30);
% 
% % Set axis limits
% xlim(x_limits);
% ylim(y_limits);
% zlim(z_limits);
% 
% % Add colorbar for time
% cb = colorbar;
% colormap(custom_cmap);
% caxis([min_delivery_time, max_delivery_time]);
% ylabel(cb, 'Time to Deliver (s)', 'FontWeight', 'bold');
% cbpos = get(cb, 'Position');
% % Move the colorbar slightly to the right
% set(cb, 'Position', [cbpos(1)+0.01 cbpos(2) cbpos(3) cbpos(4)]);
% 
% % Add legend
% h1 = plot3(NaN, NaN, NaN, 'go', 'MarkerSize', 10, 'LineWidth', 2);
% h2 = plot3(NaN, NaN, NaN, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
% h3 = plot3(NaN, NaN, NaN, 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
% 
% legend([h1, h2, h3], {'Start Point', 'Waypoint', 'Wind Vectors'}, ...
%     'Location', 'northeast');
% 
% % Box on for cleaner appearance
% box on;
% 
% %% Create data storage folders and save results
% % Create a base folder for data storage
% base_folder = 'simulation_data';
% if ~exist(base_folder, 'dir')
%     mkdir(base_folder);
% end
% 
% % Get current date in yyyyMMdd format
% current_date = datestr(now, 'yyyymmdd');
% 
% % Create date folder
% date_folder = fullfile(base_folder, current_date);
% if ~exist(date_folder, 'dir')
%     mkdir(date_folder);
% end
% 
% % Get current time for filename
% timestamp = datestr(now, 'HHMM');
% filename = sprintf('thermal_timing_comparison_%s.mat', timestamp);
% full_path = fullfile(date_folder, filename);
% 
% % Save data
% save(full_path, 'all_results', 'delivery_times', 'wind_field', 'initial_state', 'aircraft', 'waypoint');
% fprintf('\nData saved to: %s\n', full_path);
% 
% % Save figure
% fig_filename = sprintf('thermal_timing_comparison_%s.png', timestamp);
% fig_path = fullfile(date_folder, fig_filename);
% 
% % Use export_fig if available, otherwise use saveas
% try
%     export_fig(fig_path, '-m5', '-nocrop');
%     fprintf('3D flight path figure exported to: %s\n', fig_path);
% catch
%     warning('Could not export figure using export_fig, using standard print instead');
%     try
%         print(fig_path, '-dpng', '-r300');
%         fprintf('Exported using standard print function\n');
%     catch
%         warning('Failed to export figure');
%     end
% end
% 
% %% Generate summary statistics table
% summary_table = table;
% summary_table.TimeConstraint = delivery_times';
% summary_table.ActualTime = zeros(num_tests, 1);
% summary_table.MaxAltitude = zeros(num_tests, 1);
% summary_table.FinalDistance = zeros(num_tests, 1);
% summary_table.BatteryRemaining = zeros(num_tests, 1);
% summary_table.CompletedSteps = zeros(num_tests, 1);
% summary_table.ReachedWaypoint = false(num_tests, 1);
% 
% for test_idx = 1:num_tests
%     path = all_results{test_idx};
%     summary_table.ActualTime(test_idx) = path.steps * 2;
%     summary_table.MaxAltitude(test_idx) = max(path.states.z);
%     summary_table.FinalDistance(test_idx) = sqrt((waypoint.x - path.states.x(end))^2 + ...
%                                                 (waypoint.y - path.states.y(end))^2);
%     summary_table.BatteryRemaining(test_idx) = path.final_battery_percentage;
%     summary_table.CompletedSteps(test_idx) = path.steps;
% 
%     % Consider waypoint reached if within a certain distance (e.g., 50m)
%     summary_table.ReachedWaypoint(test_idx) = summary_table.FinalDistance(test_idx) < 50;
% end
% 
% % Display summary table
% disp('Summary of results:');
% disp(summary_table);
% 
% %% Create a multi-metrics plot with improved formatting
% metrics_fig = figure('Units', 'centimeters', 'Position', [0 0 17 13]);
% 
% % Set tight margins between subplots
% tight_margin = 0.12; % Normalized margin
% subplot_width = (1 - 3*tight_margin)/2;
% subplot_height = (1 - 3*tight_margin)/2;
% 
% % Plot 1: Max altitude vs time constraint (Top Left)
% ax1 = subplot(2, 2, 1);
% plot(summary_table.TimeConstraint, summary_table.MaxAltitude, 'o-', 'LineWidth', 2, 'Color', [0.2 0.5 0.7]);
% xlabel('Time Constraint (s)', 'FontWeight', 'bold');
% ylabel('Max Altitude (m)', 'FontWeight', 'bold');
% title('Maximum Altitude vs Time Constraint', 'FontWeight', 'bold');
% grid on;
% set(ax1, 'Position', [tight_margin, 1-tight_margin-subplot_height, subplot_width, subplot_height]);
% 
% % Plot 2: Battery remaining vs time constraint (Top Right)
% ax2 = subplot(2, 2, 2);
% plot(summary_table.TimeConstraint, summary_table.BatteryRemaining, 'o-', 'LineWidth', 2, 'Color', [0.8 0.4 0.2]);
% xlabel('Time Constraint (s)', 'FontWeight', 'bold');
% ylabel('Battery Remaining (%)', 'FontWeight', 'bold');
% title('Battery Remaining vs Time Constraint', 'FontWeight', 'bold');
% grid on;
% set(ax2, 'Position', [2*tight_margin+subplot_width, 1-tight_margin-subplot_height, subplot_width, subplot_height]);
% 
% % Plot 3: Actual flight time vs time constraint (Bottom Left)
% ax3 = subplot(2, 2, 3);
% plot(summary_table.TimeConstraint, summary_table.ActualTime, 'o-', 'LineWidth', 2, 'Color', [0.2 0.5 0.7]);
% hold on;
% plot(summary_table.TimeConstraint, summary_table.TimeConstraint, 'r--', 'LineWidth', 1);
% xlabel('Time Constraint (s)', 'FontWeight', 'bold');
% ylabel('Actual Flight Time (s)', 'FontWeight', 'bold');
% title('Actual vs Constrained Flight Time', 'FontWeight', 'bold');
% legend('Actual Time', 'Constraint', 'Location', 'northwest', 'FontSize', 10);
% grid on;
% set(ax3, 'Position', [tight_margin, tight_margin, subplot_width, subplot_height]);
% 
% % Plot 4: Final distance to waypoint vs time constraint (Bottom Right)
% ax4 = subplot(2, 2, 4);
% plot(summary_table.TimeConstraint, summary_table.FinalDistance, 'o-', 'LineWidth', 2, 'Color', [0.8 0.4 0.2]);
% xlabel('Time Constraint (s)', 'FontWeight', 'bold');
% ylabel('Final Distance (m)', 'FontWeight', 'bold');
% title('Distance to Waypoint vs Time Constraint', 'FontWeight', 'bold');
% grid on;
% set(ax4, 'Position', [2*tight_margin+subplot_width, tight_margin, subplot_width, subplot_height]);
% 
% % Save the metrics figure
% metrics_fig_filename = sprintf('thermal_timing_metrics_%s.png', timestamp);
% metrics_fig_path = fullfile(date_folder, metrics_fig_filename);
% 
% % Use export_fig if available, otherwise use saveas
% try
%     export_fig(metrics_fig_path, '-m5', '-nocrop');
%     fprintf('Metrics figure exported to: %s\n', metrics_fig_path);
% catch
%     warning('Could not export metrics figure using export_fig, using standard print instead');
%     try
%         print(metrics_fig_path, '-dpng', '-r300');
%         fprintf('Exported using standard print function\n');
%     catch
%         warning('Failed to export metrics figure');
%     end
% end
% 
% fprintf('\nSimulation complete. %d tests run with varying delivery time constraints.\n', num_tests);


% test_thermal_timing_comparison.m
% This script systematically tests the effect of varying delivery time 
% constraints on path planning in thermal conditions

clear all; close all; clc;

% Check if init_mission.m accepts parameters
try
    test_mission = init_mission(250);
    if test_mission.time_to_deliver ~= 250
        error('init_mission.m does not properly handle input parameters. Please update it according to the instructions.');
    end
    fprintf('init_mission.m successfully accepts time parameters.\n\n');
catch e
    error(['Error: init_mission.m needs to be modified to accept a time parameter.\n', ...
           'Please update it to have the following structure:\n', ...
           'function mission = init_mission(varargin)\n', ...
           '    mission.time_to_deliver = 300;\n', ...
           '    mission.alt_floor = 150;\n', ...
           '    mission.alt_floor_tol = 50;\n', ...
           '    if nargin > 0 && ~isempty(varargin{1})\n', ...
           '        mission.time_to_deliver = varargin{1};\n', ...
           '    end\n', ...
           'end\n']);
end

% Load configurations
aircraft = init_aircraft();
battery = init_battery();
waypoint = init_waypoint();

%% Create wind field
% Define grid to match flight space
x_limits = [-100, 2000];
y_limits = [-500, 500];
z_limits = [0, 800];
grid_spacing = 20;

% Create wind field structure
wind_field = init_wind_field(x_limits, y_limits, z_limits, grid_spacing);

% Add thermal pattern - using 'simple_sys' for consistent thermal placement
wind_field = create_thermals(wind_field, 'simple_sys');

%% Initial conditions
initial_state.x = 0;
initial_state.y = 0;
initial_state.z = 300;
initial_state.V = 15;        % m/s
initial_state.V_rel = 15;    % m/s
initial_state.gamma = 0;     % rad
initial_state.phi = 0;       % rad
initial_state.psi = 0;       % rad

%% Set up test parameters
% Define range of delivery times to test
min_delivery_time = 100;  % seconds
max_delivery_time = 400;  % seconds
num_tests = 30;

delivery_times = linspace(min_delivery_time, max_delivery_time, num_tests);
num_steps = 200;  % Maximum steps per simulation

% Storage for all test results
all_results = cell(num_tests, 1);

% Create spring colormap for time visualization
cmap = spring(num_tests);

%% Run all tests
for test_idx = 1:num_tests
    fprintf('\n\n==================================================\n');
    fprintf('STARTING TEST %d of %d: Time to deliver = %.1f seconds\n', ...
        test_idx, num_tests, delivery_times(test_idx));
    fprintf('==================================================\n\n');
    
    % Run simulation with current time parameter
    complete_path = generate_multi_step_path(initial_state, aircraft, wind_field, num_steps, battery, delivery_times(test_idx));
    
    % Store results
    all_results{test_idx} = complete_path;
    all_results{test_idx}.time_to_deliver = delivery_times(test_idx);
    
    % Calculate summary statistics for this run
    total_flight_time = complete_path.steps * 2; % seconds
    max_altitude = max(complete_path.states.z);
    final_distance = sqrt((waypoint.x - complete_path.states.x(end))^2 + ...
                         (waypoint.y - complete_path.states.y(end))^2);
    
    % Calculate final battery percentage
    battery_capacity = battery.capacity_Wh;
    base_power = battery.base_power;
    thrust_coefficient = battery.thrust_coefficient;
    
    % Simple battery calculation to estimate final percentage
    total_energy_used = 0;
    for s = 1:length(complete_path.states.thrust)
        thrust = complete_path.states.thrust(s);
        power = base_power + thrust_coefficient * thrust;
        energy = power * (2/3600); % 2 seconds converted to hours
        total_energy_used = total_energy_used + energy;
    end
    
    remaining_battery = battery_capacity - total_energy_used;
    battery_percentage = (remaining_battery / battery_capacity) * 100;
    
    all_results{test_idx}.final_battery = remaining_battery;
    all_results{test_idx}.final_battery_percentage = battery_percentage;
    
    fprintf('\nTest %d Summary:\n', test_idx);
    fprintf('  Time constraint: %.1f seconds\n', delivery_times(test_idx));
    fprintf('  Actual flight time: %.1f seconds\n', total_flight_time);
    fprintf('  Max altitude reached: %.1f meters\n', max_altitude);
    fprintf('  Final distance to waypoint: %.1f meters\n', final_distance);
    fprintf('  Remaining battery: %.2f%% (%.2f Wh out of %.2f Wh)\n', battery_percentage, remaining_battery, battery_capacity);
    fprintf('  Steps completed: %d\n', complete_path.steps);
end

%% Set up global plotting settings
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultTextInterpreter', 'tex');
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');

%% Create visualization of all paths
% Create a figure for the 3D paths with thesis-appropriate dimensions
fig_handle = figure;
set(fig_handle, 'Units', 'centimeters');
set(fig_handle, 'Position', [0 0 17 13]); % Width and height in centimeters

% Plot the wind field first (light gray)
hold on;
downsample = 3;  % Adjust for density of wind vectors
vertical_threshold = 0.2;  % Only show significant vertical winds
significant_vertical = abs(wind_field.U) >= vertical_threshold;

% Create a combined mask that both downsamples and checks for significance
mask = false(size(wind_field.X));
mask(1:downsample:end, 1:downsample:end, 1:downsample:end) = true;
mask = mask & significant_vertical;

% Plot wind vectors
quiver3(wind_field.X(mask), wind_field.Y(mask), wind_field.Z(mask), ...
        wind_field.Wx(mask), wind_field.Wy(mask), wind_field.U(mask), ...
        'Color', [0.8 0.8 0.8], 'LineWidth', 0.8, 'MaxHeadSize', 0.3);

% Plot waypoint marker
plot3(waypoint.x, waypoint.y, waypoint.z, 'r*', 'MarkerSize', 10, 'LineWidth', 2);

% Plot all paths with color coding for time to deliver
for test_idx = 1:num_tests
    path = all_results{test_idx};
    
    % Get the color for this path
    path_color = cmap(test_idx, :);
    
    % Plot this path
    plot3(path.states.x, path.states.y, path.states.z, ...
        'Color', path_color, 'LineWidth', 1.5);
    
    % Mark the start point for the first path only
    if test_idx == 1
        plot3(initial_state.x, initial_state.y, initial_state.z, ...
            'go', 'MarkerSize', 10, 'LineWidth', 2);
    end
    
    % Mark the end point of each path with a dot
    plot3(path.states.x(end), path.states.y(end), path.states.z(end), ...
        'o', 'Color', path_color, 'MarkerSize', 8, 'LineWidth', 1.5);
end

% Format plot
xlabel('East (m)', 'FontWeight', 'bold');
ylabel('North (m)', 'FontWeight', 'bold');
zlabel('Altitude (m)', 'FontWeight', 'bold');
title('Effect of Delivery Time Constraints on Path Planning', 'FontWeight', 'bold');

% Set position of axes in the figure for consistency
set(gca, 'Units', 'centimeters');
set(gca, 'Position', [2 1.5 11 10]);  % [left bottom width height]

% Set view and grid
grid off;
view(30, 30);

% Set axis limits
xlim(x_limits);
ylim(y_limits);
zlim(z_limits);

% Use actual proportions instead of axis equal
% Calculate axis ratios based on actual distances
x_range = x_limits(2) - x_limits(1);
y_range = y_limits(2) - y_limits(1);
z_range = z_limits(2) - z_limits(1);
pbaspect([x_range, y_range, z_range]);

% Add colorbar for time with discrete ticks for each test
cb = colorbar;
colormap(cmap);
caxis([min_delivery_time - (delivery_times(2)-delivery_times(1))/2, ...
      max_delivery_time + (delivery_times(2)-delivery_times(1))/2]);
% Set ticks to match the actual test points
set(cb, 'Ticks', delivery_times);
ylabel(cb, 'Time to Deliver (s)', 'FontWeight', 'bold');
cbpos = get(cb, 'Position');
% Move the colorbar slightly to the right
set(cb, 'Position', [cbpos(1)+0.01 cbpos(2) cbpos(3) cbpos(4)]);

% Add legend
h1 = plot3(NaN, NaN, NaN, 'go', 'MarkerSize', 10, 'LineWidth', 2);
h2 = plot3(NaN, NaN, NaN, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
h3 = plot3(NaN, NaN, NaN, 'Color', [0.8 0.8 0.8], 'LineWidth', 1);

legend([h1, h2, h3], {'Start Point', 'Waypoint', 'Wind Vectors'}, ...
    'Location', 'northeast');

% Box on for cleaner appearance
box on;

%% Create data storage folders and save results
% Create a base folder for data storage
base_folder = 'simulation_data';
if ~exist(base_folder, 'dir')
    mkdir(base_folder);
end

% Get current date in yyyyMMdd format
current_date = datestr(now, 'yyyymmdd');

% Create date folder
date_folder = fullfile(base_folder, current_date);
if ~exist(date_folder, 'dir')
    mkdir(date_folder);
end

% Get current time for filename
timestamp = datestr(now, 'HHMM');
filename = sprintf('thermal_timing_comparison_%s.mat', timestamp);
full_path = fullfile(date_folder, filename);

% Save data
save(full_path, 'all_results', 'delivery_times', 'wind_field', 'initial_state', 'aircraft', 'waypoint');
fprintf('\nData saved to: %s\n', full_path);

% Save figure
fig_filename = sprintf('thermal_timing_comparison_%s.png', timestamp);
fig_path = fullfile(date_folder, fig_filename);

% Use export_fig if available, otherwise use saveas
try
    export_fig(fig_path, '-m5', '-nocrop');
    fprintf('3D flight path figure exported to: %s\n', fig_path);
catch
    warning('Could not export figure using export_fig, using standard print instead');
    try
        print(fig_path, '-dpng', '-r300');
        fprintf('Exported using standard print function\n');
    catch
        warning('Failed to export figure');
    end
end

%% Generate summary statistics table
summary_table = table;
summary_table.TimeConstraint = delivery_times';
summary_table.ActualTime = zeros(num_tests, 1);
summary_table.MaxAltitude = zeros(num_tests, 1);
summary_table.FinalDistance = zeros(num_tests, 1);
summary_table.BatteryRemaining = zeros(num_tests, 1);
summary_table.CompletedSteps = zeros(num_tests, 1);
summary_table.ReachedWaypoint = false(num_tests, 1);

for test_idx = 1:num_tests
    path = all_results{test_idx};
    summary_table.ActualTime(test_idx) = path.steps * 2;
    summary_table.MaxAltitude(test_idx) = max(path.states.z);
    summary_table.FinalDistance(test_idx) = sqrt((waypoint.x - path.states.x(end))^2 + ...
                                                (waypoint.y - path.states.y(end))^2);
    summary_table.BatteryRemaining(test_idx) = path.final_battery_percentage;
    summary_table.CompletedSteps(test_idx) = path.steps;
    
    % Consider waypoint reached if within a certain distance (e.g., 50m)
    summary_table.ReachedWaypoint(test_idx) = summary_table.FinalDistance(test_idx) < 50;
end

% Display summary table
disp('Summary of results:');
disp(summary_table);

%% Create a multi-metrics plot with improved formatting
metrics_fig = figure('Units', 'centimeters', 'Position', [0 0 17 13]);

% Set tight margins between subplots
tight_margin = 0.12; % Normalized margin
subplot_width = (1 - 3*tight_margin)/2;
subplot_height = (1 - 3*tight_margin)/2;

% Plot 1: Max altitude vs time constraint (Top Left)
ax1 = subplot(2, 2, 1);
plot(summary_table.TimeConstraint, summary_table.MaxAltitude, 'o-', 'LineWidth', 2, 'Color', [0.2 0.5 0.7]);
xlabel('Time Constraint (s)', 'FontWeight', 'bold');
ylabel('Max Altitude (m)', 'FontWeight', 'bold');
title('Maximum Altitude vs Time Constraint', 'FontWeight', 'bold');
grid on;
set(ax1, 'Position', [tight_margin, 1-tight_margin-subplot_height, subplot_width, subplot_height]);

% Plot 2: Battery remaining vs time constraint (Top Right)
ax2 = subplot(2, 2, 2);
plot(summary_table.TimeConstraint, summary_table.BatteryRemaining, 'o-', 'LineWidth', 2, 'Color', [0.8 0.4 0.2]);
xlabel('Time Constraint (s)', 'FontWeight', 'bold');
ylabel('Battery Remaining (%)', 'FontWeight', 'bold');
title('Battery Remaining vs Time Constraint', 'FontWeight', 'bold');
grid on;
set(ax2, 'Position', [2*tight_margin+subplot_width, 1-tight_margin-subplot_height, subplot_width, subplot_height]);

% Plot 3: Actual flight time vs time constraint (Bottom Left)
ax3 = subplot(2, 2, 3);
plot(summary_table.TimeConstraint, summary_table.ActualTime, 'o-', 'LineWidth', 2, 'Color', [0.2 0.5 0.7]);
hold on;
plot(summary_table.TimeConstraint, summary_table.TimeConstraint, 'r--', 'LineWidth', 1);
xlabel('Time Constraint (s)', 'FontWeight', 'bold');
ylabel('Actual Flight Time (s)', 'FontWeight', 'bold');
title('Actual vs Constrained Flight Time', 'FontWeight', 'bold');
legend('Actual Time', 'Constraint', 'Location', 'northwest', 'FontSize', 10);
grid on;
set(ax3, 'Position', [tight_margin, tight_margin, subplot_width, subplot_height]);

% Plot 4: Final distance to waypoint vs time constraint (Bottom Right)
ax4 = subplot(2, 2, 4);
plot(summary_table.TimeConstraint, summary_table.FinalDistance, 'o-', 'LineWidth', 2, 'Color', [0.8 0.4 0.2]);
xlabel('Time Constraint (s)', 'FontWeight', 'bold');
ylabel('Final Distance (m)', 'FontWeight', 'bold');
title('Distance to Waypoint vs Time Constraint', 'FontWeight', 'bold');
grid on;
set(ax4, 'Position', [2*tight_margin+subplot_width, tight_margin, subplot_width, subplot_height]);

% Save the metrics figure
metrics_fig_filename = sprintf('thermal_timing_metrics_%s.png', timestamp);
metrics_fig_path = fullfile(date_folder, metrics_fig_filename);

% Use export_fig if available, otherwise use saveas
try
    export_fig(metrics_fig_path, '-m5', '-nocrop');
    fprintf('Metrics figure exported to: %s\n', metrics_fig_path);
catch
    warning('Could not export metrics figure using export_fig, using standard print instead');
    try
        print(metrics_fig_path, '-dpng', '-r300');
        fprintf('Exported using standard print function\n');
    catch
        warning('Failed to export metrics figure');
    end
end

fprintf('\nSimulation complete. %d tests run with varying delivery time constraints.\n', num_tests);