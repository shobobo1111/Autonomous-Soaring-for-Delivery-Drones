function fig_handle = plot_3d_flight_path(complete_path, wind_field, initial_state, aircraft, varargin)
    % Parse optional parameters
    p = inputParser;
    addOptional(p, 'load_file', '', @ischar);  % Filename to load
    addParameter(p, 'save_data', true, @islogical);  % Whether to save data
    parse(p, varargin{:});
    
    load_file = p.Results.load_file;
    save_data = p.Results.save_data;
    
    % Create base folder for data storage
    base_folder = 'simulation_data';
    if ~exist(base_folder, 'dir')
        mkdir(base_folder);
    end
    
    % Get current date in yyyyMMdd format using datetime
    current_time = datetime('now');
    current_date = char(string(current_time, 'yyyyMMdd'));
    
    % Create date folder
    date_folder = fullfile(base_folder, current_date);
    if ~exist(date_folder, 'dir')
        mkdir(date_folder);
    end
    
    % If load_file is provided, load data from that file
    if ~isempty(load_file)
        try
            % Check if full path or just filename
            [~, ~, ext] = fileparts(load_file);
            if isempty(ext)
                load_file = [load_file, '.mat']; % Add extension if missing
            end
            
            if ~contains(load_file, filesep)
                % If no path separator, assume it's in the current date folder
                load_file = fullfile(date_folder, load_file);
            end
            
            data = load(load_file);
            complete_path = data.complete_path;
            wind_field = data.wind_field;
            initial_state = data.initial_state;
            aircraft = data.aircraft;
            fprintf('Data loaded from: %s\n', load_file);
        catch err
            error('Error loading data from %s: %s', load_file, err.message);
        end
    else
        % Check if we need to save new data
        if save_data
            % Get current time
            % Create datetime object
            current_time = datetime('now');
            
            % Format it to string
            time_str = char(string(current_time, 'HH:mm'));
            
            % Find existing files in today's folder
            files = dir(fullfile(date_folder, '*_flight_path_data_*.mat'));
            if isempty(files)
                next_num = 1;
            else
                % Extract numbers from filenames
                nums = zeros(length(files), 1);
                for i = 1:length(files)
                    parts = split(files(i).name, '_');
                    nums(i) = str2double(parts{1});
                end
                next_num = max(nums) + 1;
            end
            
            % Create filename
            filename = sprintf('%d_flight_path_data_%s.mat', next_num, time_str);
            full_path = fullfile(date_folder, filename);
            
            % Save data
            save(full_path, 'complete_path', 'wind_field', 'initial_state', 'aircraft');
            fprintf('Data saved to: %s\n', full_path);
            
            % Store the base filename for later use in plot export
            base_filename = sprintf('%d_flight_path_%s', next_num, time_str);
        end
    end

    % Now we start plotting fr!
    % we create a 3D visualisation of the flight path with wind field
    %
    % Inputs:
    %   complete_path - Structure containing full path and state information
    %   wind_field    - Structure containing wind field data
    %   initial_state - Structure with initial aircraft state
    %   aircraft      - Structure with aircraft parameters
    %
    % Output:
    %   fig_handle    - Handle to the created figure
    
    % Global plotting settings
    set(0, 'DefaultAxesFontSize', 12);
    set(0, 'DefaultFigureColor', 'w');
    set(0, 'DefaultTextInterpreter', 'tex');
    set(0, 'DefaultAxesFontName', 'Times New Roman');
    set(0, 'DefaultTextFontName', 'Times New Roman');
    
    % Get waypoint information
    waypoint = init_waypoint();
    
    % Create figure with thesis-appropriate dimensions
    fig_handle = figure;
    set(fig_handle, 'Units', 'centimeters');
    set(fig_handle, 'Position', [0 0 17 13]); % Width and height in centimeters
    
    % Plot full wind field
    quiver3(wind_field.X, wind_field.Y, wind_field.Z, ...
        wind_field.Wx, wind_field.Wy, wind_field.U, ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 0.01);
    hold on
    
    % Define custom colormap using spectral colors
    custom_colors = {
        '#9e0142', '#d53e4f', '#f46d43', '#fdae61', '#fee08b', ...
        '#e6f598', '#abdda4', '#66c2a5', '#3288bd', '#5e4fa2'
    };
    
    % Convert hex to RGB
    custom_rgb = zeros(length(custom_colors), 3);
    for i = 1:length(custom_colors)
        hex = custom_colors{i}(2:end);
        custom_rgb(i, 1) = hex2dec(hex(1:2))/255; % Red
        custom_rgb(i, 2) = hex2dec(hex(3:4))/255; % Green
        custom_rgb(i, 3) = hex2dec(hex(5:6))/255; % Blue
    end
    
    % Interpolate to create a 100-step colormap
    custom_cmap = interp1(linspace(0, 1, size(custom_rgb, 1)), custom_rgb, linspace(0, 1, 100));
    
    % Get velocities and min/max values
    velocities = complete_path.states.V;
    vel_min = aircraft.min_speed;
    vel_max = aircraft.max_speed;
    
    % Direct mapping from velocity to colormap index
    for i = 1:length(velocities)-1
        % Map velocity directly to color index (0-100 range)
        color_idx = max(1, min(100, round((velocities(i) - vel_min) / (vel_max - vel_min) * 99) + 1));
        
        plot3([complete_path.states.x(i), complete_path.states.x(i+1)], ...
              [complete_path.states.y(i), complete_path.states.y(i+1)], ...
              [complete_path.states.z(i), complete_path.states.z(i+1)], ...
              'Color', custom_cmap(color_idx,:), 'LineWidth', 2);
    end
    
    % Plot waypoint and start point
    plot3(waypoint.x, waypoint.y, waypoint.z, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
    plot3(initial_state.x, initial_state.y, initial_state.z, 'go', 'MarkerSize', 10, 'LineWidth', 2);
    
    % Set colormap and color limits
    colormap(custom_cmap);
    clim([vel_min, vel_max]);
    
    % Add colorbar with simplified tick definition
    c = colorbar;
    ylabel(c, 'Velocity (m/s)', 'FontWeight', 'bold');
    cbpos = get(c, 'Position');
    % Move the colorbar 0.5 cm to the right
    set(c, 'Position', [cbpos(1)+0.01 cbpos(2) cbpos(3) cbpos(4)]);
    
    % Create ticks at every 5 m/s
    tick_range = vel_min:5:vel_max;
    set(c, 'Ticks', tick_range);
    
    % Create tick labels with special labels for critical speeds
    tick_labels = cell(size(tick_range));
    for i = 1:length(tick_range)
        tick_val = tick_range(i);
        if abs(tick_val - aircraft.min_speed) < 0.1
            tick_labels{i} = 'V_{min}';
        elseif abs(tick_val - aircraft.cruise_speed) < 0.1
            tick_labels{i} = 'V_{cruise}';
        elseif abs(tick_val - aircraft.max_speed) < 0.1
            tick_labels{i} = 'V_{max}';
        else
            tick_labels{i} = sprintf('%.0f', tick_val);
        end
    end
    set(c, 'TickLabels', tick_labels);
    
    % Format axes
    xlabel('East (m)', 'FontWeight', 'bold');
    ylabel('North (m)', 'FontWeight', 'bold');
    zlabel('HAGL (m)', 'FontWeight', 'bold');
    % title('Multi-Step Optimised Flight Path', 'FontWeight', 'bold');
    
    % Set position of axes in the figure for consistency
    set(gca, 'Units', 'centimeters');
    set(gca, 'Position', [2 1.5 11 10]);
    
    % Set view and grid
    grid on;
    view(-45, 30);
    axis equal;
    
    % % Add legend
    % h_waypoint = plot3(NaN, NaN, NaN, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
    % h_start = plot3(NaN, NaN, NaN, 'go', 'MarkerSize', 10, 'LineWidth', 2);
    % h_wind = plot3(NaN, NaN, NaN, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    % legend([h_waypoint, h_start, h_wind], {'Waypoint', 'Start Position', 'Wind Field'}, ...
    %     'Location', 'northeast', 'FontSize', 10);
    
    % Box off for cleaner appearance
    box off;

%     timestamp = datestr(now, 'yyyymmdd_HHMMSS');
% 
%     % For 3D flight path
%     filename = ['flight_path_3d_', timestamp, '.png'];
% 
%     % Export with high resolution
%     export_fig (filename, '-m5', '-nocrop')
%     fprintf('Figure exported as: %s\n', filename);
% end

%     % When exporting figure at the end:
%     if save_data && ~isempty(load_file)
%         % Extract base filename from loaded file
%         [~, load_name, ~] = fileparts(load_file);
%         parts = split(load_name, '_');
%         base_filename = sprintf('%s_%s_%s', parts{1}, parts{4}, parts{5});
%     end
% 
%     % Export figure using the same numbering as the data file
%     if exist('base_filename', 'var')
%         try
%             export_filename = fullfile(date_folder, [base_filename, '.png']);
%             export_fig (export_filename, '-m5', '-nocrop');
%             fprintf('Figure exported to: %s\n', export_filename);
%         catch
%             warning('Could not export figure using export_fig');
%         end
%     end
% end

% When exporting figure from loaded data:
if ~isempty(load_file)
    % Extract base number from loaded file
    [~, load_name, ~] = fileparts(load_file);
    parts = split(load_name, '_');
    original_num = parts{1};
    
    % Find existing plot files with this base number
    plot_pattern = fullfile(date_folder, [original_num, '_*_flight_path.png']);
    existing_plots = dir(plot_pattern);
    
    if isempty(existing_plots)
        sub_num = 1;
    else
        % Extract sub-numbers and find the highest
        sub_nums = zeros(length(existing_plots), 1);
        for i = 1:length(existing_plots)
            name_parts = split(existing_plots(i).name, '_');
            if length(name_parts) >= 2
                sub_nums(i) = str2double(name_parts{2});
            end
        end
        sub_num = max(sub_nums) + 1;
    end
    
    % Create filename with original number and sub-number
    export_filename = fullfile(date_folder, [original_num, '_', num2str(sub_num), '_flight_path.png']);
else
    % Normal export for new data
    export_filename = fullfile(date_folder, [base_filename, '.png']);
end

% Export the figure
try
    export_fig(export_filename, '-m5', '-nocrop');
    fprintf('Figure exported to: %s\n', export_filename);
catch
    warning('Could not export figure using export_fig');
end