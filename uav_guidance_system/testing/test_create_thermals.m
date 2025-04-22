function test_create_thermals()
    % Test visualization of different thermal types
    
    % Define grid parameters with more space around thermals
    x_limits = [-50 150];  % Expanded limits
    y_limits = [-50 150];
    z_limits = [0 200];    % Higher z limit for better thermal visualization
    grid_spacing = 10;     % Increased spacing for clearer vectors
    
    % Define thermal types to test
    thermal_types = {'simple', 'complex', 'drifting', 'twisting', ...
                    'thermal_with_sink', 'sink_drift_complex', 'simple_sys'};
    
    % Create figure with subplots
    figure('Position', [50 50 1600 900]);  % Larger figure
    
    % Loop through each thermal type
    for i = 1:length(thermal_types)
        % Create subplot with spacing between plots
        subplot(2, 4, i)
        
        % Initialize fresh wind field for each thermal
        wind_field = init_wind_field(x_limits, y_limits, z_limits, grid_spacing);
        
        % Create thermal of specific type
        wind_field = create_thermals(wind_field, thermal_types{i});
        
        % Plot 3D vector field with improved visibility
        quiver3(wind_field.X, wind_field.Y, wind_field.Z, ...
                wind_field.Wx, wind_field.Wy, wind_field.U, ...
                'Color', [0.2 0.2 0.8], ...    % Darker blue color
                'LineWidth', 1.5, ...          % Thicker lines
                'AutoScale', 'off', ...        % Disable automatic scaling
                'MaxHeadSize', 0.5);           % Larger arrow heads
            
        % Add labels and title
        xlabel('X Position (m)', 'FontSize', 12);
        ylabel('Y Position (m)', 'FontSize', 12);
        zlabel('Z Position (m)', 'FontSize', 12);
        title(['Thermal Type: ' thermal_types{i}], 'FontSize', 14);
        
        % Set consistent view angle
        view(45, 30);
        
        % Add grid and set background color
        grid on;
        set(gca, 'Color', [0.95 0.95 0.95]);  % Light grey background
        set(gca, 'GridAlpha', 0.3);           % Subtle grid
        
        % Set axis limits explicitly
        xlim(x_limits);
        ylim(y_limits);
        zlim(z_limits);
    end
    
    % Add overall title
    sgtitle('Comparison of Different Thermal Types - 3D Vector Field Visualization', ...
            'FontSize', 16);
            
    % Add space between subplots
    set(gcf, 'Units', 'normalized');
    set(gcf, 'Position', [0.1 0.1 0.8 0.8]);
end