% function plot_path_flower(paths, wind_field, step_number, elapsed_time)
%     % Create output directory if it doesn't exist
%     output_dir = 'path_flower_frames';
%     if ~exist(output_dir, 'dir')
%         mkdir(output_dir);
%     end
% 
%     % Set up figure with consistent size and position
%     fig = figure(2);
%     clf % Clear previous flower
%     set(fig, 'Position', [100 100 800 600]); % Consistent figure size
% 
%     % Plot wind field
%     quiver3(wind_field.X, wind_field.Y, wind_field.Z, ...
%         wind_field.Wx, wind_field.Wy, wind_field.U, ...
%         'k', 'LineWidth', 1);
%     hold on
% 
%     waypoint = init_waypoint();
%     plot3(waypoint.x, waypoint.y, waypoint.z, 'r*', ...
%         'MarkerSize', 10, ...
%         'LineWidth', 2);
% 
%     % Plot all paths, color coded by reward
%     rewards = [paths.total_reward];
%     min_reward = min(rewards);
%     max_reward = max(rewards);
%     reward_colors = jet(length(paths));
% 
%     % Initialize arrays to store all coordinates
%     x_all = [];
%     y_all = [];
%     z_all = [];
% 
%     for i = 1:length(paths)
%         % Plot path
%         plot3(paths(i).states.x, paths(i).states.y, paths(i).states.z, ...
%             'Color', reward_colors(i,:), 'LineWidth', 2);
%         hold on
% 
%         % Plot feeler line (dashed, lighter color)
%         feeler_color = reward_colors(i,:) * 0.8 + 0.2;
%         plot3(paths(i).feeler.x, paths(i).feeler.y, paths(i).feeler.z, ...
%             '--', 'Color', feeler_color, 'LineWidth', 1);
% 
%         % Collect coordinates for axis limits
%         x_all = [x_all; paths(i).states.x(:)];
%         y_all = [y_all; paths(i).states.y(:)];
%         z_all = [z_all; paths(i).states.z(:)];
%     end
% 
%     % Format plot
%     xlabel('X (m)')
%     ylabel('Y (m)')
%     zlabel('Z (m)')
%     title(sprintf('Path Options - Step %d (t = %.1f s)', step_number, elapsed_time))
%     grid on
% 
%     % Set consistent view angle
%     view([-37.5, 30]); % Adjust these angles as needed
% 
%     % Set axis limits based on data extent plus padding
%     x_range = [min(x_all), max(x_all)];
%     y_range = [min(y_all), max(y_all)];
%     z_range = [min(z_all), max(z_all)];
% 
%     padding = 0.2; % 20% padding around data
%     x_pad = diff(x_range) * padding;
%     y_pad = diff(y_range) * padding;
%     z_pad = diff(z_range) * padding;
% 
%     xlim([x_range(1)-x_pad, x_range(2)+x_pad]);
%     ylim([y_range(1)-y_pad, y_range(2)+y_pad]);
%     zlim([z_range(1)-z_pad, z_range(2)+z_pad]);
% 
%     % Add colorbar
%     c = colorbar;
%     ylabel(c, 'Normalised Path Reward')
%     colormap(flipud(reward_colors))
% 
%     % Print path details to command window
%     fprintf('\nPath Details (Sorted by Reward):\n')
%     fprintf('Time: %.1f s\n', elapsed_time);
%     fprintf('Path\tBank°\tPitch°\tTotal Reward\tEnd Position (x,y,z)\n')
%     for i = 1:length(paths)
%         fprintf('%d\t%.1f\t%.1f\t\t%.1f\t\t(%.1f, %.1f, %.1f)\n', ...
%             i, ...
%             paths(i).bank_angle*(180/pi), ...
%             paths(i).pitch_angle*(180/pi), ...
%             paths(i).total_reward, ...
%             paths(i).states.x(end), ...
%             paths(i).states.y(end), ...
%             paths(i).states.z(end))
%     end
% 
%     % Save figure with consistent naming
%     filename = sprintf('%s/path_flower_%04d.png', output_dir, step_number);
%     exportgraphics(fig, filename, 'Resolution', 300);
% 
%     drawnow
% end



function plot_path_flower(paths, wind_field, step_number, elapsed_time)
    % Create output directory if it doesn't exist
    output_dir = 'path_flower_frames';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % Create figure with thesis-appropriate dimensions
    fig = figure(2);
    clf
    set(fig, 'Units', 'centimeters');
    set(fig, 'Position', [0 0 17 10]); % Width and height in cm

    % Plot wind field
    quiver3(wind_field.X, wind_field.Y, wind_field.Z, ...
        wind_field.Wx, wind_field.Wy, wind_field.U, ...
        'k', 'LineWidth', 1);
    hold on

    waypoint = init_waypoint();
    plot3(waypoint.x, waypoint.y, waypoint.z, 'r*', ...
        'MarkerSize', 10, 'LineWidth', 2);

    % Plot all paths, color coded by reward
    rewards = [paths.total_reward];
    min_reward = min(rewards);
    max_reward = max(rewards);
    reward_colors = jet(length(paths));

    x_all = [];
    y_all = [];
    z_all = [];

    for i = 1:length(paths)
        plot3(paths(i).states.x, paths(i).states.y, paths(i).states.z, ...
            'Color', reward_colors(i,:), 'LineWidth', 2);
        hold on

        feeler_color = reward_colors(i,:) * 0.8 + 0.2;
        plot3(paths(i).feeler.x, paths(i).feeler.y, paths(i).feeler.z, ...
            '--', 'Color', feeler_color, 'LineWidth', 1);

        x_all = [x_all; paths(i).states.x(:)];
        y_all = [y_all; paths(i).states.y(:)];
        z_all = [z_all; paths(i).states.z(:)];
    end

    % Format axes
    xlabel('East (m)', 'FontWeight', 'bold');
    ylabel('North (m)', 'FontWeight', 'bold');
    zlabel('HAGL (m)', 'FontWeight', 'bold');
    title(sprintf('Path Options - Step %d (t = %.1f s)', step_number, elapsed_time), ...
        'FontWeight', 'normal');

    % Set consistent view angle
    view(45, 10);

    % Axis limits with padding
    padding = 0.2;
    x_range = [min(x_all), max(x_all)];
    y_range = [min(y_all), max(y_all)];
    z_range = [min(z_all), max(z_all)];

    x_pad = diff(x_range) * padding;
    y_pad = diff(y_range) * padding;
    z_pad = diff(z_range) * padding;

    xlim([x_range(1)-x_pad, x_range(2)+x_pad]);
    ylim([y_range(1)-y_pad, y_range(2)+y_pad]);
    zlim([z_range(1)-z_pad, z_range(2)+z_pad]);

    % Layout adjustment with precise padding
    padding_left = 0.2;
    padding_bottom = 1;
    padding_right = 0.2;
    padding_top = 0.1;

    fig_width = 17;
    fig_height = 10;
    plot_width = fig_width - padding_left - padding_right;
    plot_height = fig_height - padding_bottom - padding_top;

    set(gca, 'Units', 'centimeters');
    set(gca, 'Position', [padding_left, padding_bottom, plot_width, plot_height]);

    % Colorbar
    c = colorbar;
    ylabel(c, 'Normalised Path Reward', 'FontWeight', 'bold');
    colormap(flipud(reward_colors));
    
    % Turn off grid if preferred
    grid off;

    % Print path details
    fprintf('\nPath Details (Sorted by Reward):\n');
    fprintf('Time: %.1f s\n', elapsed_time);
    fprintf('Path\tBank°\tPitch°\tTotal Reward\tEnd Position (x,y,z)\n');
    for i = 1:length(paths)
        fprintf('%d\t%.1f\t%.1f\t\t%.1f\t\t(%.1f, %.1f, %.1f)\n', ...
            i, ...
            paths(i).bank_angle*(180/pi), ...
            paths(i).pitch_angle*(180/pi), ...
            paths(i).total_reward, ...
            paths(i).states.x(end), ...
            paths(i).states.y(end), ...
            paths(i).states.z(end));
    end

    % Save high-res figure
    filename = sprintf('%s/path_flower_%04d.png', output_dir, step_number);
    exportgraphics(fig, filename, 'Resolution', 300);

    drawnow
end
