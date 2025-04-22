% test_path_generation.m
clear all; close all; clc;

%% Import aircraft parameters
aircraft = init_aircraft();

%% Create wind field
% Define grid to match flight space
x_limits = [-100, 300];    % Plenty of room for paths to develop
y_limits = [-100, 300];    % Symmetric around start point
z_limits = [0, 400];       % From ground to above flight level
grid_spacing = 20;         % Coarse enough to see vectors clearly

% Create wind field structure
wind_field = init_wind_field(x_limits, y_limits, z_limits, grid_spacing);

% Add desired thermal pattern
wind_field = create_thermals(wind_field, 'simple');

%% Test 1: Generate paths from level flight
% Initial conditions
current_state.x = 40;
current_state.y = 10;
current_state.z = 200;
current_state.V = 15;      % m/s
current_state.V_rel = 10;      % m/s
current_state.gamma = 0.0;   % rad
current_state.phi = 0;     % rad
current_state.psi = 0;     % rad

elapsed_time = 0;

operating_mode = 1;

paths = generate_paths(current_state, aircraft, wind_field, elapsed_time, operating_mode);

%% Plot results
figure(1)
% Plot wind field first
quiver3(wind_field.X, wind_field.Y, wind_field.Z, ...
    wind_field.Wx, wind_field.Wy, wind_field.U, ...
    'k', 'LineWidth', 1);
hold on

% Plot all paths and their feelers, color coded by reward
rewards = [paths.total_reward];
min_reward = min(rewards);
max_reward = max(rewards);
[custom_cmap, named_colors] = init_colour();
reward_colors = flipud(custom_cmap(round(linspace(1, size(custom_cmap, 1), length(paths))), :));

for i = 1:length(paths)
    % Plot main path
    plot3(paths(i).states.x, paths(i).states.y, paths(i).states.z, ...
        'Color', reward_colors(i,:), 'LineWidth', 2)
    
    % Plot feeler line (dashed, lighter color)
    if isfield(paths(i), 'feeler') && isfield(paths(i).feeler, 'x')
        feeler_color = reward_colors(i,:) * 0.8 + 0.2;  % Make color lighter
        plot3(paths(i).feeler.x, paths(i).feeler.y, paths(i).feeler.z, ...
            '--', 'Color', feeler_color, 'LineWidth', 1)
    end
end

% Format plot
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
title('Paths and Feeler Lines Colored by Total Reward')
grid on
view(3)
axis equal

% Add colorbar for rewards
c = colorbar;
ylabel(c, 'Path Reward')
colormap(custom_cmap)

% Add legend
h_path = plot3(NaN,NaN,NaN, '-k', 'LineWidth', 2);
h_feeler = plot3(NaN,NaN,NaN, '--k', 'LineWidth', 1);
legend([h_path, h_feeler], {'Main Path', 'Feeler Line'}, 'Location', 'best')

%% Print path details with rewards
fprintf('\nPath Details (Sorted by Reward):\n')
fprintf('Path\tBank°\tPitch°\tPath Reward\tFeeler Reward\tTotal Reward\tEnd Position (x,y,z)\tV_rel\tV\n')
for i = 1:length(paths)
    fprintf('%d\t%.1f\t%.1f\t%.1f\t\t%.1f\t\t%.1f\t\t(%.1f, %.1f, %.1f)\t%.1f\t%.1f\n', ...
        i, ...
        paths(i).bank_angle*(180/pi), ...
        paths(i).pitch_angle*(180/pi), ...
        paths(i).total_reward - paths(i).feeler_reward, ...
        paths(i).feeler_reward, ...
        paths(i).total_reward, ...
        paths(i).states.x(end), ...
        paths(i).states.y(end), ...
        paths(i).states.z(end), ...
        paths(i).states.V_rel(end), ...
        paths(i).states.V(end))
end

%% Additional Top-Down View with Graying Feelers
figure(2)
% Plot wind field first (top down view of vectors)
quiver(wind_field.X(:,:,1), wind_field.Y(:,:,1), ...
    wind_field.Wx(:,:,1), wind_field.Wy(:,:,1), ...
    'b', 'LineWidth', 1);
hold on

% Get feeler config
feeler_config = init_feeler_config();
num_feeler_segments = feeler_config.num_points - 1;

% Plot paths and feelers
rewards = [paths.total_reward];
min_reward = min(rewards);
max_reward = max(rewards);
[custom_cmap, named_colors] = init_colour();
reward_colors = custom_cmap(round(linspace(1, size(custom_cmap, 1), length(paths))), :);


for i = 1:length(paths)
    % Plot main path
    plot(paths(i).states.x, paths(i).states.y, ...
        'Color', reward_colors(i,:), 'LineWidth', 2)
    
    % Plot feeler line with darkening gray
    if feeler_config.enabled
        if isfield(paths(i), 'feeler') && isfield(paths(i).feeler, 'x')
            for f = 1:num_feeler_segments
                % Calculate gray shade based on decay
                gray_level = 0.7 * (1 - (feeler_config.decay_rate^f));  
                
                % Plot feeler segment
                plot([paths(i).feeler.x(f), paths(i).feeler.x(f+1)], ...
                     [paths(i).feeler.y(f), paths(i).feeler.y(f+1)], ...
                     '--', 'Color', [gray_level gray_level gray_level], ...
                     'LineWidth', 1)
            end
        end
    end
end

% Format plot
xlabel('X (m)')
ylabel('Y (m)')
title('Top-Down View with Decaying Feeler Lines')
grid on
axis equal

% Add colorbar and legend
c = colorbar;
ylabel(c, 'Path Reward')
colormap(reward_colors)
h_path = plot(NaN,NaN, '-b', 'LineWidth', 2);
h_feeler = plot(NaN,NaN, '-b', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
legend([h_path, h_feeler], {'Main Path', 'Feeler (Darker = greater Influence)'}, ...
    'Location', 'best')

%% ==== Plot Config and High-Res Export ====
% Define plot and export settings
plot_width_cm = 17;             % Width in cm
aspect_ratio = 4/3;             % Default aspect ratio (change if needed)
plot_height_cm = plot_width_cm / aspect_ratio;
dpi = 600;                     % 10x resolution if default is 300
font_size = 12;                 % Default font size
padding_cm = 1.5;               % Padding around plot for labels

% Set figure size in pixels for export
width_inch = plot_width_cm / 2.54;
height_inch = plot_height_cm / 2.54;

set(gcf, 'Units', 'inches', 'Position', [1, 1, width_inch, height_inch])
set(gca, 'FontSize', font_size)

% Export as high-resolution PNG
export_name = 'path_reward_plot';
print(gcf, export_name, '-dpng', ['-r', num2str(dpi)]);
disp(['Exported figure to ', export_name, '.png at ', num2str(dpi), ' DPI'])
